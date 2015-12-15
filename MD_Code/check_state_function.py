def check_state_function(coordinates):
    if -100.0 < coordinates[0] < 0.0 and -120.0 < coordinates[1] < 30.0 and -100.0 < coordinates[2] < 0.0 and \
       -120.0 < coordinates[3] < 30.0 and -100.0 < coordinates[4] < 0.0 and -120.0 < coordinates[5] < 30.0:
        return 0
    elif -180.0 < coordinates[0] < -55.0 and (105.0 < coordinates[1] <= 180.0 or -180.0 <= coordinates[1] < -155.0) and \
         -180.0 < coordinates[2] < -55.0 and (105.0 < coordinates[3] <= 180.0 or -180.0 <= coordinates[3] < -155.0) and \
         -180.0 < coordinates[4] < -55.0 and (105.0 < coordinates[5] <= 180.0 or -180.0 <= coordinates[5] < -155.0):
        return 1
    else:
        return -1
