# Synaptic-Self-Organization-of-Spatio-Temporal-Pattern-Selectivity


 This repository contains the codes to implement the model described in the manuscript "**Synaptic Self-Organization of Spatio-Temporal Pattern Selectivity**".



**To reproduce the results, follow the steps below:**

Extract the Folder.zip into the desired directory. You will see a folder called *Folders*. Then, locate all the MATLAB Files (.m files) into the *Folders*. 



**Fig1, Fig2, Fig3, Fig4, and Fig5:**

First execute Model_Learning.m for p1=1:1:500, and p2=1:1:3.
Then to plot figure 1, execute  Fig1.m, figure 2a: Fig_2a.m, and so on. (Note: For Fig4, change Ntrials_3=20000 to Ntrials_3=60000 )



**Fig6:**

First, execute then Model_Learning.m for  p1=1:1:500, and p2=1:1:3. second execute Model_test_Noise.m for  p1=1:1:500, and p2=1:1:3. Then execute Fig_6a.m, Fig_6b.m, and Fig_6c.m



**Fig7:**

First: execute the Model_Learning_noise.m for  p1=1:1:500 and p2=1. 
Second: execute Model_test_Noise_2.m for p1=1:1:500 and p2=1.
Third: execute  Fig_7a.m, Fig_7b.m, and Fig_7c.m



**Fig8a,b ,and d:**

First execute Model_Learning_2_em for p1=1:500 and p2=1. Then execute Fig_8ab.m and Fig_8d.m



**Fig8c:**

First, execute the Model_Learning_4_em.m for  p1=1:1:100 and p2=1. 
Second execute Fig_8c.m



**Fig9a:**

First, execute Model_Learning_NET_1.m (p1:1:50 and p2=1:4), then Fig_9a.m


**Fig9b, c:**

First, execute Model_Learning_NET_2.m (p1:1:500 and p2=4), then Fig_9b.m, Fig_9c.m 



**Fig9d:**

This figure was obtained from an extensive simulation. For every number of embedded patterns, we considered post-synaptic neurons from 2:2:30. Each point is from an average of 50 simulations.
Model_Learning_NET_1.m is the primary program to plot Fig 9d. The positions of the embedded patterns are fixed; however, results are independent of the pattern position.


**Fig10:**

First: Execute Model_Learning_NET_3.m (p1:1:50)
Second: execute Fig_10.m



**Fig12:**

First: execute the Model_Learning_f0.m for  p1=1:1:100 and p2=1:3. 
Second: execute Fig_12.m

