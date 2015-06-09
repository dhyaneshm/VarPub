import multiprocessing

def worker(lines):
    """Make a dict out of the parsed, supplied lines"""
    #result = {}
    result = 0
    #lines = '-'.join(lines)
    print str(lines)

    for line in lines.split('\n'):
        result = result + int(line)
    #for line in lines.split('\n'):
    #    k, v = parse(line)
    #    result[k] = v
    return result

if __name__ == '__main__':
    # configurable options.  different values may work better.
    numthreads = 2
    numlines = 100

    lines = open('input.txt').readlines()

    # create the process pool
    pool = multiprocessing.Pool(processes=numthreads)

    # map the list of lines into a list of result dicts
    #print str(len(lines)) + " " + str(type(lines) is list)

    #for line in xrange(0,len(lines),numlines):
    #    print lines[line:line+numlines]

    result_list = pool.map(worker, ( lines[line:line+numlines] for line in xrange(0,len(lines),numlines) )  )

    print '\n'.join(result_list)

    # reduce the result dicts into a single dict
    result = {}
    map(result.update, result_list)

