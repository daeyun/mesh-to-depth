import numpy as np
import pyassimp


def read_mesh_assimp(filename):
    scene = pyassimp.load(filename)

    vertices = []
    faces = []
    vertex_count = 0
    for m in scene.meshes:
        f = m.faces
        v = m.vertices
        faces.append(f + vertex_count)
        vertices.append(v)
        num_vertices = v.shape[0]
        vertex_count += num_vertices
    faces = np.concatenate(faces, axis=0)
    vertices = np.concatenate(vertices, axis=0)

    fv = {'v': vertices, 'f': faces}

    pyassimp.release(scene)

    return fv
