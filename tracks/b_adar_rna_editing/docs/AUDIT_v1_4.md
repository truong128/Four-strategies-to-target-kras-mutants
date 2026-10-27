# v1.4 adversarial audit

## Verified in the container

- 21 Python unit and regression tests pass.
- Syntax compilation for all Python source/tests and Bash run script passes.
- Synthetic E2E invocation and input-required behavior exercised by existing tests.
- Bowtie2 index isolation enforced for all configured transcriptome runs.
- Positional chemistry inter-nucleotide PO/PS and terminal bond semantics tested.

## Limitations and unverified areas

- External command interface compatibility with installed GATK/VEP/Boltz/Amber/OpenMM releases is not verified by these tests.
- Published RESTORE chemistry transfer does not constitute KRAS editing validation.
- Force-field fidelity, ADAR catalytic geometry and measured efficacy require independently curated input and/or empirical experiments.
- An external independent test suite and release validation on target machine remain essential.

## Reproducible checks

```bash
./run_all.sh smoke
python3 -m unittest discover -s tests -p 'test_*.py' -v
python3 -m compileall -q full scripts tests
bash -n run_all.sh
```
