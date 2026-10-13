% FILENAME:  mybatch.m

pjob=batch('measure_phase_delay4bands','Profile','SlurmBell1','Pool',100,'CaptureDiary',true);
wait(pjob);
diary(pjob);
quit;
