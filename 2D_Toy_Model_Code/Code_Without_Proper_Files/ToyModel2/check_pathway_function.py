def check_pathway_function(coordinates):
    if -0.75 <= coordinates[0] <= 0.75 and coordinates[1] < 0.25:
        return 0
    elif -0.75 <= coordinates[0] <= 0.75 and 0.25 <= coordinates[1]:
        return 1
    else:
        return -1
