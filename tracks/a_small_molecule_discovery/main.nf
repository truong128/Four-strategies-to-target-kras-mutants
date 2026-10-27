nextflow.enable.dsl=2

params.config = params.config ?: 'config.yaml'
params.library = params.library ?: null
params.training_docking_csv = params.training_docking_csv ?: null
params.receptor_pdbqt = params.receptor_pdbqt ?: null
params.center_x = params.center_x ?: null
params.center_y = params.center_y ?: null
params.center_z = params.center_z ?: null
params.size_x = params.size_x ?: 22
params.size_y = params.size_y ?: 22
params.size_z = params.size_z ?: 22
params.outdir = params.outdir ?: 'results'
params.run_admet = params.run_admet == null ? true : params.run_admet
params.receptor_pdb = params.receptor_pdb ?: 'structures/7RPZ_KRAS_G12D_protein.pdb'

process STANDARDIZE_LIBRARY {
  publishDir "${params.outdir}/library", mode: 'copy'
  cpus 4; memory '8 GB'
  input:
    path library
    path cfg
  output:
    path 'standardized.smi', emit: standardized
    path 'library_qc.json', emit: qc
  script:
  """
  python ${projectDir}/scripts/prepare_library.py --input ${library} --output standardized.smi \
    --config ${cfg} --qc library_qc.json
  """
}

process SAMPLE_LIBRARY {
  publishDir "${params.outdir}/sample", mode: 'copy'
  cpus 28; memory '16 GB'
  input:
    path library
    path cfg
  output:
    path 'training_sample.smi'
  script:
  """
  python ${projectDir}/scripts/sample_library.py --input ${library} --output training_sample.smi --config ${cfg}
  """
}

process DOCK_TRAINING {
  publishDir "${params.outdir}/docking", mode: 'copy'
  cpus 28; memory '16 GB'
  input:
    path sample
    path receptor
    path cfg
  output:
    path 'training_docking.csv'
  script:
  """
  python ${projectDir}/scripts/dock_vina_batch.py --smiles ${sample} --receptor ${receptor} \
    --output training_docking.csv --config ${cfg} \
    --center ${params.center_x} ${params.center_y} ${params.center_z} \
    --size ${params.size_x} ${params.size_y} ${params.size_z}
  """
}

process TRAIN_ML {
  publishDir "${params.outdir}/ml", mode: 'copy'
  cpus 12; memory '32 GB'
  input:
    path docking
    path cfg
  output:
    path 'model_comparison.csv', emit: comparison
    path 'model_bundle.joblib', emit: bundle
    path 'locked_test_metrics.json', emit: metrics
    path 'split_manifest.csv', emit: split
  script:
  """
  python ${projectDir}/scripts/benchmark_16_models.py --input ${docking} --config ${cfg} \
    --comparison model_comparison.csv --bundle model_bundle.joblib \
    --metrics locked_test_metrics.json --split-manifest split_manifest.csv
  """
}

process SCREEN_LIBRARY {
  publishDir "${params.outdir}/screen", mode: 'copy'
  cpus 28; memory '32 GB'
  input:
    path library
    path bundle
    path cfg
  output:
    path 'ml_hits.csv.gz'
  script:
  """
  python ${projectDir}/scripts/predict_ultralarge.py --input ${library} --bundle ${bundle} \
    --config ${cfg} --output ml_hits.csv.gz
  """
}

process REDOCK {
  publishDir "${params.outdir}/redock", mode: 'copy'
  cpus 28; memory '16 GB'
  input:
    path hits
    path receptor
    path cfg
  output:
    path 'redocked_hits.csv'
  script:
  """
  python ${projectDir}/scripts/dock_vina_batch.py --csv ${hits} --receptor ${receptor} \
    --output redocked_hits.csv --config ${cfg} \
    --center ${params.center_x} ${params.center_y} ${params.center_z} \
    --size ${params.size_x} ${params.size_y} ${params.size_z} \
    --redock
  """
}

process ADMET {
  publishDir "${params.outdir}/admet", mode: 'copy'
  cpus 8; memory '24 GB'
  input:
    path hits
    path cfg
  output:
    path 'admet_predictions.csv', emit: predictions
    path 'admet_filtered.csv', emit: filtered
  script:
  """
  python ${projectDir}/scripts/run_admet.py --input ${hits} --config ${cfg} \
    --predictions admet_predictions.csv --filtered admet_filtered.csv
  """
}

process CLUSTER {
  publishDir "${params.outdir}/cluster", mode: 'copy'
  cpus 8; memory '24 GB'
  input:
    path hits
    path cfg
  output:
    path 'diverse_hits.csv'
  script:
  """
  python ${projectDir}/scripts/cluster_hits.py --input ${hits} --config ${cfg} --output diverse_hits.csv
  """
}

process SYNTHESIS {
  publishDir "${params.outdir}/synthesis", mode: 'copy'
  cpus 4; memory '32 GB'
  input:
    path hits
    path cfg
    path boltz_done
  output:
    path 'synthesis_scores.csv'
  script:
  """
  python ${projectDir}/scripts/synthetic_feasibility.py --input ${hits} --config ${cfg} \
    --output synthesis_scores.csv
  """
}

process BOLTZ {
  publishDir "${params.outdir}/boltz", mode: 'copy'
  cpus 4; memory '32 GB'
  input:
    path hits
    path cfg
    path protein
  output:
    path 'boltz_scores.csv'
  script:
  """
  python ${projectDir}/scripts/run_boltz_batch.py --input ${hits} --protein-pdb ${protein} \
    --config ${cfg} --output-dir boltz_out
  cp boltz_out/boltz_scores.csv boltz_scores.csv
  """
}

process MD {
  publishDir "${params.outdir}/md", mode: 'copy'
  cpus 4; memory '32 GB'
  input:
    path hits
    path cfg
    path protein
    path retro_done
  output:
    path 'md_scores.csv'
  script:
  """
  python ${projectDir}/scripts/run_md_batch.py --input ${hits} --protein ${protein} --config ${cfg} \
    --output md_scores.csv --workdir md_work
  """
}

process MMPBSA {
  publishDir "${params.outdir}/mmpbsa", mode: 'copy'
  cpus 4; memory '8 GB'
  input:
    path hits
    path cfg
    path protein
  output:
    path 'mmpbsa_scores.csv'
  script:
  """
  python ${projectDir}/scripts/run_mmpbsa_batch.py --input ${hits} --protein ${protein} --config ${cfg} \
    --output mmpbsa_scores.csv --workdir mmpbsa_work
  """
}

process PREPARE_RESISTANCE_PANEL {
  publishDir "${params.outdir}/resistance", mode: 'copy'
  cpus 4; memory '16 GB'
  input:
    path cfg
    path panel_tsv
  output:
    path 'panel', emit: panel
  script:
  """
  python ${projectDir}/scripts/build_resistance_panel.py --config ${cfg} --panel ${panel_tsv} --output-dir panel
  """
}

process DOCK_RESISTANCE {
  publishDir "${params.outdir}/resistance", mode: 'copy'
  cpus 16; memory '16 GB'
  input:
    path hits
    path panel
    path cfg
  output:
    path 'resistance_scores.csv', emit: scores
    path 'resistance_long.csv', emit: long
  script:
  """
  python ${projectDir}/scripts/dock_resistance_panel.py --input ${hits} --panel ${panel} \
    --config ${cfg} --output resistance_scores.csv --long resistance_long.csv
  """
}

process FINAL_RANK {
  publishDir "${params.outdir}/final", mode: 'copy'
  cpus 2; memory '8 GB'
  input:
    path base
    path admet
    path synth
    path boltz
    path md
    path mmp
    path resistance
    path cfg
  output:
    path 'ranked_candidates.csv'
    path 'pareto_candidates.csv'
  script:
  """
  python ${projectDir}/scripts/final_rank.py --base ${base} --admet ${admet} --synthesis ${synth} \
    --boltz ${boltz} --md ${md} --mmpbsa ${mmp} --resistance ${resistance} \
    --config ${cfg} --ranked ranked_candidates.csv --pareto pareto_candidates.csv
  """
}

workflow {
  if (!params.library) error "Provide --library"
  cfg = Channel.value(file(params.config))
  rawLib = file(params.library)
  streamLib = rawLib.isDirectory() || params.library ==~ /.*\.(bz2|gz)$/
  lib = Channel.value(rawLib)

  if (streamLib) {
    sampleLib = lib
    screenLib = lib
  } else {
    STANDARDIZE_LIBRARY(lib, cfg)
    STANDARDIZE_LIBRARY.out.standardized.tap { sampleLib }.set { screenLib }
  }
  SAMPLE_LIBRARY(sampleLib, cfg)

  if (params.training_docking_csv) {
    docking = Channel.fromPath(params.training_docking_csv, checkIfExists: true)
  } else {
    if (!params.receptor_pdbqt || params.center_x == null || params.center_y == null || params.center_z == null)
      error "Provide receptor_pdbqt and docking-box center, or training_docking_csv"
    receptor = Channel.fromPath(params.receptor_pdbqt, checkIfExists: true)
    DOCK_TRAINING(SAMPLE_LIBRARY.out, receptor, cfg)
    docking = DOCK_TRAINING.out
  }

  TRAIN_ML(docking, cfg)
  SCREEN_LIBRARY(screenLib, TRAIN_ML.out.bundle, cfg)

  if (params.receptor_pdbqt && params.center_x != null && params.center_y != null && params.center_z != null) {
    receptor2 = Channel.fromPath(params.receptor_pdbqt, checkIfExists: true)
    REDOCK(SCREEN_LIBRARY.out, receptor2, cfg)
    base = REDOCK.out
  } else {
    base = SCREEN_LIBRARY.out
  }

  if (params.run_admet) {
    ADMET(base, cfg)
    CLUSTER(ADMET.out.filtered, cfg)
    admet_for_rank = ADMET.out.predictions
  } else {
    CLUSTER(base, cfg)
    admet_for_rank = base
  }

  CLUSTER.out.tap { forSynth }.tap { forBoltz }.tap { forMd }.tap { forMmp }.tap { forResist }.set { forRank }
  protein = Channel.fromPath(params.receptor_pdb, checkIfExists: true)
  panelTsv = Channel.fromPath("${projectDir}/config/pdb_panel.tsv", checkIfExists: true)
  PREPARE_RESISTANCE_PANEL(cfg, panelTsv)
  BOLTZ(forBoltz, cfg, protein)
  BOLTZ.out.tap { boltzForRetro }.set { boltzForRank }
  SYNTHESIS(forSynth, cfg, boltzForRetro)
  SYNTHESIS.out.tap { synthForMd }.set { synthForRank }
  MD(forMd, cfg, protein, synthForMd)
  MMPBSA(forMmp, cfg, protein)
  DOCK_RESISTANCE(forResist, PREPARE_RESISTANCE_PANEL.out.panel, cfg)
  FINAL_RANK(forRank, admet_for_rank, synthForRank, boltzForRank, MD.out, MMPBSA.out, DOCK_RESISTANCE.out.scores, cfg)
}
