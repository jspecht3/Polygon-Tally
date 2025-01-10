import csv

# parsing fission
ids = []
vals = []

with open('get_fission_out.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[0] == '0':
            continue
        if row[0] == 'time':
            ids = row[1:]
        if row[0] == '1':
            vals = row[1:]

fission_dict = {}

for i in range(len(ids)):
    _, ring, cell = ids[i].split('_')
    ring = int(ring)
    cell = int(cell)

    if ring not in fission_dict:
        fission_dict[ring] = {}

    fission_dict[ring][cell] = float(vals[i])

fission_vals = {}

for ring in fission_dict:
    if ring == 0:
        cell_count = 1
    if ring != 0:
        cell_count = ring * 6

    fission_vals[ring] = []

    for i in range(cell_count):
        fission_vals[ring].append(fission_dict[ring][i])

# parsing flux
ids = []
vals = []

with open('get_flux_out.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[0] == '0':
            continue
        if row[0] == 'time':
            ids = row[1:]
        if row[0] == '1':
            vals = row[1:]

flux_dict = {}

for i in range(len(ids)):
    _, ring, cell = ids[i].split('_')
    ring = int(ring)
    cell = int(cell)

    if ring not in flux_dict:
        flux_dict[ring] = {}

    flux_dict[ring][cell] = float(vals[i])

flux_vals = {}

for ring in flux_dict:
    if ring == 0:
        cell_count = 1
    if ring != 0:
        cell_count = ring * 6

    flux_vals[ring] = []

    for i in range(cell_count):
        flux_vals[ring].append(flux_dict[ring][i])


def count_pins():
    count = 0
    for ring in flux_vals:
        for cell in flux_vals[ring]:
            count += 1
    print(f"Total of {count} pins")
    return


# count_pins()
