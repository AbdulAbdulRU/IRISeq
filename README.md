# IRISeq
IRISeq—an optics-free, cost-effective platform that leverages spatial interaction mapping by indexed sequencing to profile tissues at adjustable sizes and resolutions (5–50 µm). 

<img width="2260" height="1288" alt="image" src="https://github.com/user-attachments/assets/be3e2146-6718-4af5-84fb-313f36274d6c" />


In this repository, we provide IRISeq pipeline for aligning cDNA reads, and also for beads beads interactions. Moreover, we provide two pipelines for tissue reconstructions. GPU based for fast and efficient reconstructions for high resolution and large areas arrays, and also CPU based for smaller ~50 um arrays. 

We also provide detailed code analysis for the IRISeq manuscript 

In Branch cDNA_Processing, Bead_interaction_pipeline and EasySpatial folder contain all related code to align reads for cDNAs from both tissue RNA and bead-beads interactions. 
The scripts to run these pipeline are cDNA_processing_script and beads_interaction_processing_script both in Branch_cDNA_Processing
