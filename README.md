# Project
Thesis preparation project: "Improving prediction models on metabarcoding data"

ABSTRACT:

Predictive modelling is a tool in data analysis to predict an unknown class or predict a future
event with statistical methods. Applying predictive models to metabarcoding data takes the
underlying patterns of present and absent species, and then predict different kinds of
ecological factors. This method was used in the publication of ​ Fløjgaard et al.’s (2019) to
predict environmental factors, geographical locations and different habitat types. The data
samples used was metabarcoded eDNA from soil and came from 130 different sites across
Denmark. ​ Fløjgaard et al ​ .’s (2019) study aimed to show the potential for provenance
prediction in forensic science and the importentness of eDNA in forensic ecology. While
showing good results for predicting environmental gradients, the prediction of geographic
origins and habitat types showed some difficulties.
This project improves the used prediction models, and therefore can predict the ecological
habitats and geographic locations more reliably. The results were improved by changing the
prediction techniques; the Random Forest and the K-nearest neighbor classifier were added
to the already used Quadratic discriminant analysis. Furthermore, multiple predictors were
added. To be able to trust the predictions more, the metric was changed from Accuracy to
Kappa. The publications model for the classification of the Forest habitat was implemented
to allow for comparison.
The best results for Kappa were shown by having predictors from the fungi ordination data
and the plant ordination data (model 5). The kNN method computed a Kappa of 0.920, where
the recreated model showed a Kappa of 0.889.
Model 5 was used on more unbalanced classification of the Atlantic geographical origin and
the Oak and Willow habitat. The results showed an improvement of the value of true
positives.
Ultimately, better results are shown when changing the classification method and adding
more predictor variables, especially predictor variables from different ecological species.

