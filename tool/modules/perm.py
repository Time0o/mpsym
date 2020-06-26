def cycle_form(generator):
    cycle = []
    cycles = []

    first = 1
    current = 1
    done = set()

    while True:
        done.add(current)
        cycle.append(current)

        current = generator[current - 1]

        if current == first:
            if len(cycle) > 1:
                cycles.append(cycle[:])

            if len(done) == len(generator):
                break

            cycle.clear()

            for i in range(1, len(generator) + 1):
                if i not in done:
                    first = i
                    current = i
                    break

    if not cycles:
        return '()'

    def print_cycle(cycle):
        return '(' + ','.join([str(x) for x in cycle]) + ')'

    return ''.join([print_cycle(cycle) for cycle in cycles])
