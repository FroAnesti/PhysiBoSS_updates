# PhysiBoSS_updates
All the changes made in the PhysiBoSS app

The less changes in their app the better. So, although, firstly, I was doing the normalization in the app, then I decided to do
the normalization in my program. Thus, the only thing I changed in their app is to create the producer, contruct the message
and send the message to Kafka.

pid: Acquire this simulation's pid using get_pid().
timepoint: the time in which this file is refering to.
aliveNO, apoptoticNO, NecroticNO: count the number of alive, apoptotic and necrotic cells in each timepoint.
path: where to find this simulation
