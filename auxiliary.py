import itertools
import re
import pathlib
import random

random.seed(42)
# print(list(itertools.permutations([0, 1, 2, 3, 4])))

# {{0,1,2,3,4},{0,1,2,4,3},{0,1,3,2,4},{0,1,3,4,2},
#  {0,1,4,2,3},{0,1,4,3,2},{0,2,1,3,4},{0,2,1,4,3},
#  {0,2,3,1,4},{0,2,3,4,1},{0,2,4,1,3},{0,2,4,3,1},
#  {0,3,1,2,4},{0,3,1,4,2},{0,3,2,1,4},{0,3,2,4,1},
#  {0,3,4,1,2},{0,3,4,2,1},{0,4,1,2,3},{0,4,1,3,2},
#  {0,4,2,1,3},{0,4,2,3,1},{0,4,3,1,2},{0,4,3,2,1},
#  {1,0,2,3,4},{1,0,2,4,3},{1,0,3,2,4},{1,0,3,4,2},
#  {1,0,4,2,3},{1,0,4,3,2},{1,2,0,3,4},{1,2,0,4,3},
#  {1,2,3,0,4},{1,2,3,4,0},{1,2,4,0,3},{1,2,4,3,0},
#  {1,3,0,2,4},{1,3,0,4,2},{1,3,2,0,4},{1,3,2,4,0},
#  {1,3,4,0,2},{1,3,4,2,0},{1,4,0,2,3},{1,4,0,3,2},
#  {1,4,2,0,3},{1,4,2,3,0},{1,4,3,0,2},{1,4,3,2,0},
#  {2,0,1,3,4},{2,0,1,4,3},{2,0,3,1,4},{2,0,3,4,1},
#  {2,0,4,1,3},{2,0,4,3,1},{2,1,0,3,4},{2,1,0,4,3},
#  {2,1,3,0,4},{2,1,3,4,0},{2,1,4,0,3},{2,1,4,3,0},
#  {2,3,0,1,4},{2,3,0,4,1},{2,3,1,0,4},{2,3,1,4,0},
#  {2,3,4,0,1},{2,3,4,1,0},{2,4,0,1,3},{2,4,0,3,1},
#  {2,4,1,0,3},{2,4,1,3,0},{2,4,3,0,1},{2,4,3,1,0},
#  {3,0,1,2,4},{3,0,1,4,2},{3,0,2,1,4},{3,0,2,4,1},
#  {3,0,4,1,2},{3,0,4,2,1},{3,1,0,2,4},{3,1,0,4,2},
#  {3,1,2,0,4},{3,1,2,4,0},{3,1,4,0,2},{3,1,4,2,0},
#  {3,2,0,1,4},{3,2,0,4,1},{3,2,1,0,4},{3,2,1,4,0},
#  {3,2,4,0,1},{3,2,4,1,0},{3,4,0,1,2},{3,4,0,2,1},
#  {3,4,1,0,2},{3,4,1,2,0},{3,4,2,0,1},{3,4,2,1,0},
#  {4,0,1,2,3},{4,0,1,3,2},{4,0,2,1,3},{4,0,2,3,1},
#  {4,0,3,1,2},{4,0,3,2,1},{4,1,0,2,3},{4,1,0,3,2},
#  {4,1,2,0,3},{4,1,2,3,0},{4,1,3,0,2},{4,1,3,2,0},
#  {4,2,0,1,3},{4,2,0,3,1},{4,2,1,0,3},{4,2,1,3,0},
#  {4,2,3,0,1},{4,2,3,1,0},{4,3,0,1,2},{4,3,0,2,1},
#  {4,3,1,0,2},{4,3,1,2,0},{4,3,2,0,1},{4,3,2,1,0}}


permutations = list(itertools.permutations([0, 1, 2, 3, 4]))  # permutations of [x(i), y(i), z(i), w(i), u(i)]; we need to do pairwise comparisons to create the lookup tabble
# pairwise comparisons for the lookup tablbe: 
# x(i) > y(i), x(i) > z(i), x(i) > w(i), x(i) > u(i)
# y(i) > z(i), y(i) > w(i), y(i) > u(i)
# z(i) > w(i), z(i) > u(i)
# w(i) > u(i)


# from 4D: c = (x0 > y0)*32 + (x0 > z0)*16 + (y0 > z0)*8 + (x0 > w0)*4 + (y0 > w0)*2 + (z0 > w0) <- for clarity: x0 == x(i), "i" being the permutation index
# in 5D:
# (x(i) > y(i))*512 + (x(i) > z(i))*256 + (x(i) > w(i))*128 + (x(i) > u(i))*64 + (y(i) > z(i))*32 + (y(i) > w(i))*16 + (y(i) > u(i))*8 + (z(i) > w(i))*4 + (z(i) > u(i))*2 + (w(i) > u(i))*1

simplex5: list[tuple[int, ...]] = [(0, 0, 0, 0, 0) for _ in range(1024)]

# comparisons translate to:
# simplex5entryindex = [permutationsentry[0] > permutationsentry[1], permutationsentry[0] > permutationsentry[2], permutationsentry[0] > permutationsentry[3], permutationsentry[0] > permutationsentry[4],
#  permutationsentry[1] > permutationsentry[2], permutationsentry[1] > permutationsentry[3], permutationsentry[1] > permutationsentry[4],
#  permutationsentry[2] > permutationsentry[3], permutationsentry[2] > permutationsentry[4],
#  permutationsentry[3] > permutationsentry[4]]

for permutationsindex, permutationsentry in enumerate(permutations):
    simplex5entryindextable = [permutationsentry[0] > permutationsentry[1], permutationsentry[0] > permutationsentry[2], permutationsentry[0] > permutationsentry[3], permutationsentry[0] > permutationsentry[4],
                          permutationsentry[1] > permutationsentry[2], permutationsentry[1] > permutationsentry[3], permutationsentry[1] > permutationsentry[4],
                          permutationsentry[2] > permutationsentry[3], permutationsentry[2] > permutationsentry[4],
                          permutationsentry[3] > permutationsentry[4]]
    simplex5entryindex = sum([simplex5entryindextable[len(simplex5entryindextable) - 1 - index] * 2**index for index in range(len(simplex5entryindextable) - 1, 0 - 1, -1)])
    simplex5[simplex5entryindex] = permutationsentry

counter = 0
for elem in simplex5:
    if elem != (0, 0, 0, 0, 0):
        counter += 1

assert len(permutations) == counter

simplex5string = re.sub(r"[\[|(]", "{", str(simplex5))
simplex5string = re.sub(r"[\]|)]", "}", simplex5string)
simplex5string = re.sub(r" ", "", simplex5string)
simplex5string = re.sub("},", "}, ", simplex5string)
simplex5stringlist = re.split(r" ", simplex5string)

index = 0
for counter in range(len(simplex5stringlist)):
    if not (counter + 1) % 8:
        # print(-(counter + 1))
        index += 1
        simplex5stringlist.insert(index + counter, "\n")

simplex5string = "".join(simplex5stringlist)
simplex5string = "const unsigned char SIMPLEX5[][5] = {\n" + simplex5string[1:-1] + ";\n\n"

perms = random.choices(range(0, 255), k=1024)

permsstring = str(perms)
permsstring = re.sub(r"[\[]", "{", permsstring)
permsstring = re.sub(r"[\]]", "}", permsstring)

permsstringlist = permsstring.split(",")

index = 0
for counter in range(len(perms)):
    if not (counter + 1) % 32:
        # print(-(counter + 1))
        index += 1
        permsstringlist.insert(index + counter, "\n")

permsstring = ",".join(permsstringlist)
permsstring = re.sub(" ", "", permsstring)
permsstring = re.sub("\n,", "\n", permsstring)
permsstring = "const unsigned char PERM5[] = {\n" + permsstring[1:-2] + ";\n\n"


gradperms = list(itertools.product([1, -1], repeat=4))
numperms = len(gradperms)

grad5 = [gradperms[index % numperms][0: int(index // numperms)] + (0,) + gradperms[index % numperms][int(index // numperms):] for index in range(5 * numperms)]
assert len(grad5) == 80


grad5string = re.sub(r"[\[|(]", "{", str(grad5))
grad5string = re.sub(r"[\]|)]", "} ", grad5string)
grad5string = re.sub(" ", "", grad5string)
grad5string = re.sub("},", "}, ", grad5string)
grad5stringlist = re.split(" ", grad5string)

index = 0
for counter in range(len(grad5stringlist)):
    if not (counter + 1) % 8:
        # print(-(counter + 1))
        index += 1
        grad5stringlist.insert(index + counter, "\n")

grad5string = "".join(grad5stringlist)

grad5string = "const float GRAD5[][5] = {\n" + grad5string[1:-1] + ";"

savepath = pathlib.Path(__file__).parent

with open(f"{savepath}\\_noise5d.h", "w") as file:
    inpt = simplex5string + grad5string
    file.write(inpt)
