# Intraparietal stimulation disrupts negative distractor effects in human multi-alternative decision-making

Kohl, C.<sup>1</sup>, Wong, M. X. M. <sup>1</sup>, Rushworth, M. F. S. <sup>2</sup> & Chau, B. K. H. <sup>1, 3</sup>  
   
<sup><sup>1</sup> Department of Rehabilitation Sciences, The Hong Kong Polytechnic University, Hong Kong  
<sup>2</sup> Department of Experimental Psychology, University of Oxford, Oxford OX1 3UD, UK  
<sup>3</sup> University Research Facility in Behavioral and Systems Neuroscience, The Hong Kong Polytechnic University, Hong Kong</sup>  

***
This repository contains behavioural data as well as code to replicate all main behavioural findings associated with the manuscript “Intraparietal stimulation disrupts negative distractor effects in human multi-alternative decision-making”. 

***
## Code
This repository contains three .m files which run all main behavioural analyses, using the data provided in the ‘Data’ directory.
* **GLM1.m**
    *	This script fits GLM1 to all NonTMS trials   
    *	GLM1:	β0 + β1 z(HV-LV) + β2 z(HV+LV) + β3 z(D-HV) + β4 z(HV-LV) z(D-HV) + ε   
    *   One-sample t-tests on the resulting beta weights show a negative (HV-LV)(D-HV) effect on accuracy   
    *	This script plots Figure 2c   
    
*	**GLM2.m**
    *	This script fits GLM2 to all NonTMS trials, after applying a HV-LV median split
    *	GLM2:	Step 1, β0 + β1 z(HV-LV) + β2 z(HV+LV) + ε1 - Step 2, β3 + β4 z(D-HV) + ε2
    *   One-sample t-tests on the resulting beta weights show a positive distractor effect on hard trials and a negative distractor effect on easy trials
    *	This script plots Figure 2d
    
*	**GLM3.m**
    *	This script fits GLM3 separately to each condition ( MIP/MT x TMS/NonTMS x Contralateral D/Ipsilateral D)
    *	GLM3:	β0 + β1 z(HV-LV) + β2 z(HV+LV) + β3 z(D-HV) + ε
    *   The resulting beta weights of each regressor are entered in a Session x Stimulation x Distractor Location ANOVA 
    *	This script plots Figure 3a-c
    

***
## Data
The ‘Data’ directory contains one .mat file per participant (01-31). Each file contains a structure called ‘data’ with three fields:
*	data.<span>MIP</span> contains data collected in the MIP session (270 x 18, rows = trials)
*	data.<span>MT</span> contains data collected in the MT session (270 x 18, rows = trials)
*	data.<span>Key</span> contains a cell with 18 strings, labelling each column in <span>data</span>.MT/MIP       

    <sub>1: **Trial_Nr**	Trial number *(1-270)*   
    2: **HV_value**	Value of the high-value option *(value = magnitude x probability)*   
    3: **LV_value**	Value of the low-value option *(value = magnitude x probability)*   
    4: **D_value**	Value of the distractor *(value = magnitude x probability)*   
    5: **HV_magnitude**	Magnitude of the high-value option *(1-6)*      
    6: **LV_magnitude**	Magnitude of the low-value option *(1-6)*    
    7: **D_magnitude**	Magnitude of the distractor *(1-6)*    
    8: **HV_magnitude**	Probability of reward associated with the high-value option *(12.5-87.5)*   
    9: **LV_magnitude**	Probability of reward associated with the low-value option *(12.5-87.5)*   
    10: **D_magnitude**	Probability of reward associated with the distractor *(12.5-87.5)*   
    11: **HV_position**	Position of the high-value option *(1 = top left, 2 = top right, 3 = bottom left, 4 = bottom right)*   
    12: **LV_position**	Position of the low-value option *(1 = top left, 2 = top right, 3 = bottom left, 4 = bottom right)*   
    13: **D_position**	Position of the distractor *(1 = top left, 2 = top right, 3 = bottom left, 4 = bottom right)*   
    14: **TMS**	TMS applied *(1 = TMS, 0 = NonTMS)*   
    15: **Decision** Choice made *(1 = HV chosen, 2 = LV chosen, nan = D/empty quadrant chosen)*   
    16: **RT**	Reaction time in ms *(max. 1500)*   
    17: **Reward**	Reward achieved   
    18: **Accuracy**  Choice accuracy *(1 = HV chosen, 0 = LV chosen, na n= D/empty quadrant chosen)* </sub>   

***


Further information, code, and data may be available upon request. 
Please refer to the manuscript or contact: kohl.carmen.1@gmail.com or boltonchau@gmail.com
