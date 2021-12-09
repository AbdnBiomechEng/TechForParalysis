# -*- coding: utf-8 -*-

"""
Test with pytrigno to connect to localhost and get emg data
"""

import time
import sys
sys.path.append("./pytrigno/")
try:
    import pytrigno
except ImportError:
    import sys
    sys.path.insert(0, '..')
    import pytrigno

# Start the EMG
dev = pytrigno.TrignoEMG(channel_range=(0, 7), samples_per_read=1,host='127.0.0.1')
# test multi-channels
dev.set_channel_range((0, 7)) # connect to the 8 channels of the emg
dev.start()
print("Started Trigno")

try:
	while True:
		print('\nAsking for emg data')
		data = dev.read() # get data from trigno
		# we get a numpy.darray (8, 270) # print(type(data), data.shape)
		datalist = data.tolist() # instead of serialize/json
		time.sleep(0.1)
		print(datalist)
		
except KeyboardInterrupt:
	print("\nAsked to stop trigno")
	dev.stop()
	print("Stopped trigno")
