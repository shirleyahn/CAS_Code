def check_pathway_function(coordinates):
    if -0.5 <= coordinates[0] <= 0.5 and coordinates[1] <= 0.25:
        return 0
    elif -0.5 <= coordinates[0] <= 0.5 and 0.5 <= coordinates[1]:
        return 1
    else:
        return -1
