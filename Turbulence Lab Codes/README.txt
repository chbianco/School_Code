

xwire_calibrate_acquire:   acquire data for xwire calibration, saves raw calibration data suitable for....

xwire_calibrate_process:   standalone code to process evaluate and save calibration data for xwires

OR

xwire_compute_calibration: function that can be called from another script to compute xwire calibration

once you have the calibration: 

xwire_takedata:  record hotwire data at a number of velocities, save the raw voltages in a .mat file; plot the power spectrum and the raw voltages on top of the calibration data.  Saves the voltages AND the calibration data

xwire_processdata:  take the raw voltages and convert them to u and v. Plot the spectra and the distribution of u' and v'.  Saves the voltages and velocities and the calibration data