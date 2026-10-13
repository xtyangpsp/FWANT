% FILENAME:  mybatch.m

pjob=batch('kern_synthetic4kernel','Profile','SlurmBell1','Pool',30,'CaptureDiary',true);
wait(pjob);
diary(pjob);
quit;
