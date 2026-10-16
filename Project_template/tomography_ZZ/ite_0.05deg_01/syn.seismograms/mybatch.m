% FILENAME:  mybatch.m

pjob=batch('snap2sac_spher','Profile','SlurmBell1','Pool',30,'CaptureDiary',true);
wait(pjob);
diary(pjob);
quit;
