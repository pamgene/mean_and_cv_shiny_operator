# Mean and CV Shiny_operator

##### Description

The `mean_and_cv_shiny_operator` is a shiny app used for creating CV plots, and creating a fit of the variation as function of intensity to the Two Component Error Model.

##### Usage

Input projection|.
---|---
`row`   | represents the variables (e.g. ID)
`col`   | represents the category (e.g. barcode, Sample Name)
`y-axis`| measurement value
`color`	| color (e.g. barcode, Sample Name)

Output relations|.
---|---
`Operator view`        | view of the Shiny application

##### Details

The operator is a shiny app for creating an CV, and SNR and an SD plot. It can also display the error model fit inside the plot. Next to this, there are inputs to set the quantile and the x and y axis limits. 

##### See Also

[mean_operator](https://github.com/tercen/mean_operator)
[mean_sd_operator](https://github.com/tercen/mean_sd_operator)
