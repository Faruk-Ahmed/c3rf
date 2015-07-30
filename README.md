###Code for optimizing expected IoU using Candidate Constrained CRFs.  
  
Project page: https://filebox.ece.vt.edu/~faruk/c3rf/main.html  
  
* ./intseg  
  * code: Code for performing binary segmentation.   
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;          Run demo_intseg_c3rf.  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;           - Setting radius = 1 is CRF-EIoEU+enum.  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;           - Cross-validation over intermediate values is C^3RF-EIoEU+enum.  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;         To replicate the plots we report, (edit and) run plot_results.m.  
  * data: Graph structures, node and edge potentials.  
* ./pascalseg  
  * code: Code for performing semantic segmentation on PASCAL VOC 2012 val.  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;         Run demo_pascalseg_c3rf.  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;           - Setting radius = 1 is CRF-EIoEU+enum.  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;           - Cross-validation over intermediate values is C^3RF-EIoEU+enum.  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;         To replicate the plots we report, (edit and) run plot_results.m.  
  * data: Graph structures, node and edge potentials, precomputed diverse   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;          solutions as created by Premachandran et al.   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;          You will need to download the VOCdevkit folder, and edit   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;          demo_pascalseg_c3rf.m in ../code/ to include the path for the   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;          folder.  

Mex files might need recompiling.  
Feel free to get in touch at faruk -at- vt -dot- edu if you have questions or comments or if you just want to talk to me.
