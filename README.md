# Gene Expression Classification: Influenza A Effect on pDCs
--------------------------------------------------------------

Overview
--------------------------------------------------------------
This project analyzes gene expression data to study the effects of Influenza A 
infection on plasmacytoid dendritic cells (pDCs). Using the dataset 
‘GSE30550_series_matrix.txt’, it applies preprocessing, ANOVA-based feature 
selection, dimensionality reduction, and classification modeling to identify 
patterns distinguishing infected vs. control samples.

The pipeline includes:
1. Data loading and cleaning
2. Preprocessing and normalization
3. ANOVA F-test for feature selection (top 200 genes)
4. Model training and evaluation with cross-validation
5. PCA visualization
6. Heatmap of top genes

--------------------------------------------------------------
Dependencies
--------------------------------------------------------------
This script requires the following Python libraries:

- numpy
- pandas
- seaborn
- matplotlib
- scikit-learn
- xgboost
- lightgbm
- (optional) gseapy

To install all dependencies:
    pip install numpy pandas seaborn matplotlib scikit-learn xgboost lightgbm gseapy

--------------------------------------------------------------
Input Data — Downloading the Dataset
--------------------------------------------------------------
Dataset: **GSE30550 — Influenza A virus infection of human plasmacytoid dendritic cells (pDCs)**

You can obtain the data in two ways:

**Option 1: Manual Download**
1. Go to the NCBI GEO dataset page:  
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30550  
2. Scroll to the “Download family” section.  
3. Download the file named **GSE30550_series_matrix.txt.gz**
4. Extract the file (e.g., using WinRAR, 7-Zip, or `gunzip` on macOS/Linux):
   gunzip GSE30550_series_matrix.txt.gz
5. Place the extracted file (`GSE30550_series_matrix.txt`) in the same directory 
   as the Python script (`Influenza.py`).

**Option 2: Command Line Download (macOS/Linux)**
If you have `wget` installed:
   wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30550/matrix/GSE30550_series_matrix.txt.gz
   gunzip GSE30550_series_matrix.txt.gz

Then confirm the file exists with:
   ls GSE30550_series_matrix.txt

--------------------------------------------------------------
How to Run
--------------------------------------------------------------
1. Make sure the dataset is available in the working directory.
2. Open a terminal in this folder and run:

    python Influenza.py

3. The script will:
   - Load and preprocess the expression matrix
   - Perform ANOVA feature selection (top 200 genes)
   - Train and evaluate multiple classifiers using cross-validation
   - Display model performance metrics and confusion matrices
   - Visualize PCA and gene expression heatmaps

--------------------------------------------------------------
Models Evaluated
--------------------------------------------------------------
The script compares the following classifiers:

- Logistic Regression
- Support Vector Machine (RBF Kernel)
- Random Forest
- Gradient Boosting
- K-Nearest Neighbors (KNN)
- XGBoost
- LightGBM

Each model is evaluated using 5-fold stratified cross-validation.

Metrics reported:
- Accuracy
- F1 Score (weighted)
- ROC-AUC (one-vs-rest)

--------------------------------------------------------------
Output & Visualization
--------------------------------------------------------------
The following visual outputs are generated:

- Confusion matrix for each model
- PCA scatter plot of top 200 ANOVA-selected genes
- Heatmap showing expression of selected genes
- Cross-validation summary table printed in the console

If GSEApy is installed, pathway enrichment analysis is also attempted.

--------------------------------------------------------------
Notes
--------------------------------------------------------------
- Some classifiers (like LightGBM or XGBoost) may require extra memory; 
  if you encounter issues, reduce the number of folds or disable parallel jobs.
- The dataset labels are inferred heuristically based on sample names 
  ("hour", "post", or "flu" → Infected).

--------------------------------------------------------------
Example Results
--------------------------------------------------------------
After running, you will see outputs similar to:

    Selected top 200 genes using ANOVA F-test
    Evaluating Random Forest...
    CV Accuracy: 0.94 | F1: 0.93 | ROC-AUC: 0.95
    Test Accuracy: 0.92 | F1: 0.91 | ROC-AUC: 0.94
    Classification Report:
                precision  recall  f1-score  support
      Control       0.90     0.93     0.91
      Infected      0.94     0.92     0.93


