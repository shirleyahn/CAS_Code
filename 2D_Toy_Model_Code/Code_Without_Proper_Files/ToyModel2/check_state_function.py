def check_state_function(coordinates):
    if -1.25 <= coordinates[0] <= -0.75 and -0.05 <= coordinates[1] <= 0.25:
        return 0
    elif 0.75 <= coordinates[0] <= 1.25 and -0.05 <= coordinates[1] <= 0.25:
        return 1
    else:
        return -1
