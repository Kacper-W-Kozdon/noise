import itertools
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
print(len(simplex5))

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

print(simplex5[0:10])

# for i in range(3 - 1, 0 - 1, -1):
#     print(i)
