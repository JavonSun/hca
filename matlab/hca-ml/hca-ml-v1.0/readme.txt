Function HcaMl (HcaMl.m) implements the maximum likelihood method for heritable component analysis, which takes family pedigree data as genetic input. The existence of related individuals in the data is required to apply this method. 

Reference papers:
An Effective Method to Identify Heritable Components from Multivariate Phenotypes
Jiangwen Sun, Henry R Kranzler and Jinbo Bi
PloS One, 10 (12), 2015

Quadratic Optimization to Identify Highly Heritable Quantitative Traits from Complex Phenotypic Features
Jiangwen Sun, Jinbo Bi and Henry R Kranzler
In the Proceedings of the 19th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, 2013


Function get_families (get_families.m) is a utility function, which extracts families from given pedigree file and organizes family members in parent - progeny order. This function is designed for preparing input for function PedHQTI and function kinship_matrix (see below).


Function kinship_matrix (kinship_matrix.m) is also a utility function. This function creates the kinship matrix, which is one of the inputs of function PedHQTI.


Function match (match.m) is a utility function that locates elements in first input vector from the second input vector.


For more usage info (such as detailed input and output specifications) of each function, please refer to the header of the specific file that contains the function implementation.


An analysis example with both script and data is provided in the example folder


The manuscript submitted to PlosOne is included in folder plosone-manuscript
