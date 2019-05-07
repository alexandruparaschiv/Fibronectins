#!/usr/bin/env python3


# breadth-first search algo for clustering the fibrils

import queue
import pdb

test_matrix = [ [0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0] ]
test_matrix = [ [0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0] ]


def bfs(matrix):

    visited = [ False for _ in range(len(matrix)) ]
    fibrils = []
    q = queue.Queue()

    for i in range(len(matrix)-1):
        if visited[i]:
            continue
        else:
            fibrils.append([i])
            q.put(i)
        while not q.empty():
            to_visit = q.get()
            if not visited[to_visit]:
                visited[to_visit] = True
                for j in range(len(matrix)):
                    if matrix[to_visit][j]:
                        if not visited[j]:
                            q.put(j)
                            fibrils[-1].append(j);visited[j]=True
    return fibrils

print(bfs(test_matrix))
