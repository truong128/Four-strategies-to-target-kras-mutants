# 20 classifiers

The training step fits 20 models on the docked sample and writes `model_comparison.csv`.

Sixteen names follow the classifier set of the author's earlier tabular benchmark (BRFSS, `get_classifiers`):

1. `logistic_regression`
2. `multinom_reg`
3. `svm_rbf`
4. `svm_linear`
5. `decision_tree`
6. `random_forest`
7. `gradient_boosting` (histogram gradient boosting, the CPU model in that repository)
8. `extra_trees`
9. `adaboost`
10. `bag_tree`
11. `bagging`
12. `knn`
13. `lda`
14. `qda` (eigen solver with shrinkage, so it can fit 2048-bit fingerprints)
15. `naive_bayes`
16. `mlp`

Four tabular neural nets:

17. `ft_transformer` - feature-token transformer, 3 layers, 64-d tokens, bits grouped by 32
18. `ft_transformer_lite` - 2 layers, 32-d tokens, bits grouped by 64
19. `realmlp` - scaled 3-layer SELU MLP
20. `tabular_resnet` - residual MLP, width 128, 3 blocks

Model selection uses validation average precision. The locked test is scored once.

Only `production_eligible` models in `config.yaml` can become the scorer for the multi-billion library. That list is the models whose prediction cost stays inside the planned cycle: logistic regression, multinomial logistic, gradient boosting, naive Bayes, LDA, a single tree, the sklearn MLP, RealMLP, and Tabular ResNet. SVMs, k-nearest neighbors, QDA, large tree ensembles, and both FT-Transformers are still trained and appear in the comparison table. They are not used to scan the full library.
