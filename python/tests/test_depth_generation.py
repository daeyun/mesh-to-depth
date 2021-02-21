import os
import json
import numpy as np
import trimesh

import mesh_to_depth as m2d

resources_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', 'resources'))


def test_mesh2depth_simple():
    with open(os.path.join(resources_dir, 'airplane/transforms_train.json'), 'r') as f:
        cameras = json.load(f)

    near = cameras['camera_nearClip']
    far = cameras['camera_farClip']
    fov = cameras['camera_angle_x']

    mesh = trimesh.load(os.path.join(resources_dir, 'airplane/models/model_normalized.obj'), force='mesh')

    v = (np.array(mesh.vertices) + cameras['object_movement']).astype(np.float32)
    f = np.array(mesh.faces).astype(np.uint32)

    params = []
    for i in range(3):
        params.append({
            'cam_pos': cameras['frames'][i]['camera_pos'],
            'cam_lookat': cameras['frames'][i]['lookAt_pos'],
            'cam_up': [0, 1, 0],
            'x_fov': fov,
            'near': near,
            'far': far,
            'height': 256,
            'width': 256,
            'is_depth': True,
        })

    depth_maps = m2d.mesh2depth(v, f, params, empty_pixel_value=np.nan)
    for p in params:
        p['is_depth'] = False
    ray_displacement_maps = m2d.mesh2depth(v, f, params, empty_pixel_value=np.nan)

    assert len(depth_maps) == 3
    assert len(ray_displacement_maps) == 3

    for i in range(len(depth_maps)):
        assert np.isnan(depth_maps[i]).sum() == np.isnan(ray_displacement_maps[i]).sum()
        valid_pixel_ratio = np.isfinite(depth_maps[i]).sum() / np.prod(depth_maps[i].shape)
        assert 0.1 < valid_pixel_ratio < 0.3

        valid_values = depth_maps[i][np.isfinite(depth_maps[i])]
        assert valid_values.min() > 3
        assert valid_values.max() < 5

        valid_values = ray_displacement_maps[i][np.isfinite(ray_displacement_maps[i])]
        assert valid_values.min() > 3
        assert valid_values.max() < 5

        valid_values = (depth_maps[i] - ray_displacement_maps[i])[np.isfinite(depth_maps[i])]
        assert (valid_values < 1e-7).all()
        assert (valid_values > -0.5).all()
