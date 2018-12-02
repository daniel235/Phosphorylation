# breastCancer
protein data on different types of breast cancer


Topic: Bioinformatics and Computational Biology
Cases of the disease are often considered a single outcome, and assumed to share a common etiology.
However, evidence indicates that many human diseases arise and evolve through a range of
heterogeneous molecular pathologic processes. That’s why one drug can improve the symptoms of a
patient while making the other patient with the same disease worse. For this reason, even for one disease,
there might be multiple subtypes (subcategores). For example, there are five main molecular subtypes of
breast cancer. For example, Luminal which tend to grow slowly and have the best prognosis and Basal
which is more common to happen with a mutation in a gene called BRCA1 and it is more common
between younger women and have worse prognosis.
In this problem, different molecular data of 16 breast cancer patients have been measures. 8 of these
patients are in the Basal category and 8 of them are in the Luminal category.
Datasets: Available at http://bit.ly/UTRGVHackR
- Protein expression: Proteins are the smallest functional unit in the body. Each protein has a
function. Sometimes proteins are working in a group to do a job (i.e. pathways). In this dataset,
the intensity of expression of proteins in different patient are measured and reported.
- Protein-Interaction Network: This dataset tells you what proteins are interacting with each other.
Physical interaction of proteins might reflect the functional similarities between them
- Protein Phosphorylation: Phosphorylation is a modification that happens to a protein. It is
basically a phosphate compound that is attached to a specific parts of proteins. Phosphorylation
has been shown to effect lots of signaling and essential pathways in the body. Most of the drugs
that have been developed are preventing unnecessary phosphorylation of proteins. In this dataset,
the intensity of how much each protein locations has been phosphorylated for each patient has
been measured and reported. (Each protein might have multiple different location/site to be
phosphorylated)
- Kinase-Substrate Associations: Kinases are the protein who are responsible to attach the
phosphate compound to the protein. But each kinase has specific property that will be able to
phosphorylate specific group of proteins. In this data set, which protein (i.e. kinase)
phosphorylate specific part of which protein (i.e. substrate) has been reported.
Challenge:
Given all these data, can you find interesting pattern in the data? The following are some examples:
- Can you find specific pattern/feature that distinguish basal from luminal?
- Can you find any relationship among these data? For example can you say proteins that are
interacting are more likely to expressed highly together? Can you say proteins that are interacting
are more likely to phosphorylate together? Can you say the sites that are on the same protein are
more likely to phosphorylate together? How can you show these relationships? How can you
assess the significance of these relationships?
- Can you cluster the kinases or the substrates of kinases in a meaningful way to find some pattern
in their expression or phosphorylation?
- Can you identify subset of phosphorylated site or proteins that distinguish basal from luminal?
How can you assess the significance of this relationship?
