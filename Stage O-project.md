
### The curse of Dimensionality: A Basic Introduction <a id="The curse of Dimensionality-a basic introduction"></a>

The world of data science and machine learning is centered around the use of data and algorithm. In recent time there has been surge in accumulation of  data termed "Big Data". Inefficiency in handling big data could lead to inefficient result.The attributes and variables used to build a data model can enhance or reduced the accuracy of the model.
A dataset with large numbers of attributes is referred to as "High dimensional Data". such could be a curse or a breakthrough. This essay focuses on the curse of dimensionality.

### The curse of Dimensionality <a id="the curse of dimensionality"><a>
Historically, "the term curse of dimensionality" was first used by Richard E. Bellman back in the 1960s. Bellman stated that "the number of samples needed to estimate an arbitrary function with a given level of accuracy grows exponentially with respect to the number of the input variables(ie dimensionality) of the function. An increase in the data 
dimension will cause the data space to grow leading to challenges in analysis of the data and machine models less accurate.
in the world of data science as it relates to buildings models to aid cancer treatment and research, many variables and attributes are utilized during data collection. Most of the variables also has sub variables.This attributes varies from clinical attributes(patient demographic, medical history, vital sign etc),genomic and molecular attributes 
(Gene expression, mutations, protein expression levels etc), Imaging data, reatment related variables and a host number of other variables.


### Effect of the curse of Dimensionality <a id="Effect of the curse of diemnsionality"></a>
There are several effect that is associated with the curse of dimensionality on the machine model one of such is that there will be complexity in data analysis less accuracy in the prediction of the machine model. In cancer research there is data generated technology like genomics and proteomics, each sample will generate thousands of attributes with
such, the efficiency of the statistical  model and visualization becomes reduced.
Again there is a problem  of overfitting risk where the model can work efficiently on trained data, but poorly on test data.
Data Sparsity is another effect of the curse of dimensionality, with the surge in data volume it might be had to detect meaningful patterns because the data point are spread thinly across high dimensionaly space.

### Combating the Effect of Dimensionality <a id="Combating the Effect of Dimensionality" ></a>
The first step in combating the effect of dimensionality is dimension reduction. A process called PCA (Principal Component Analysis) can be used in Gene Expression Profiling where a data scientist interprete and visualize (using t.SNT) in lower dimension high dimensional data while preserving the local structure of the data and is able to differentiate 
between cancerous tissue and non-cancerous tissue based on the gene expressions patterns.
Also the process of feature selection, slicing a subset of the attributes from the original data can be used. One of such method is LASSO (Least Absolute Shrinkage and Operator). In cancer research it focuses on most informatives genetic makers.

### Conclusion<a id="conclusion"></a>
The Application of Data Science and Machine Learning model to Cancer Research is a vigorous and data packed process. The goal of a model is to have a higher accuracy to predict that which it is model for, Hence the topic of dimensionality of attributes is a MUST as a better understanding will eventually lead to High accuracy of the model.
