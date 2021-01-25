/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <iostream>
#include <fstream>

NORI_NAMESPACE_BEGIN

void OctTree::build(const BoundingBox3f& bouding_box, const MeshArray &meshes)
{
    m_log_file = fopen("test.log", "w");
    m_mesh_array = meshes;
    // 将mesh TriangeIndex 取出来
    size_t tri_cnt = 0;
    for(auto mesh_itr=meshes.begin(); mesh_itr != meshes.end(); mesh_itr++)
    {
        tri_cnt += (*mesh_itr)->getTriangleCount();
    }

    TriArray triangles;
    triangles.reserve(tri_cnt);
    for (size_t mesh_idx=0; mesh_idx<meshes.size(); mesh_idx++)
    {
        Mesh* mesh = meshes[mesh_idx];
        for (size_t tri_idx=0; tri_idx<mesh->getTriangleCount(); tri_idx++)
        {
            TriangleIndex triangle(mesh_idx, tri_idx);
            triangles.push_back(triangle);
        }
    }
    root = make_tree(bouding_box, triangles, meshes);
}

OctTreeNodePtr OctTree::make_tree(const BoundingBox3f& bounding_box, TriArray triangles, const MeshArray &meshes)
{
    if(triangles.size() <= 0)
    {
        return nullptr;
    }
    
    OctTreeNodePtr node_ptr = std::make_shared<OctTreeNode>();

    // 如果triangles的数量小于10 则创建leaf node
    size_t leaf_cnts = 100;
    if(triangles.size() <= leaf_cnts)
    {
        node_ptr->bounding_box = bounding_box;
        node_ptr->m_tri_arr = triangles;
        return node_ptr;
    }
    // 生成8个子节点和8个和子节点相交的triangle list，递归调用 make_tree
    Point3f min_pt = bounding_box.min;
    Point3f max_pt = bounding_box.max;
    BoundingBox3f sub_boxes[8];

    BoundingBox3f tmp_box;

    Point3f half_pt = 0.5 * (min_pt + max_pt);

    Point3f tmp_min = min_pt;
    Point3f tmp_max = half_pt;
    tmp_box.reset();
    tmp_box.expandBy(tmp_min);
    tmp_box.expandBy(tmp_max);
    sub_boxes[0] = tmp_box;

    tmp_min = Point3f(half_pt.x(), min_pt.y(), min_pt.z());
    tmp_max = Point3f(max_pt.x(), half_pt.y(), half_pt.z());
    tmp_box.reset();
    tmp_box.expandBy(tmp_min);
    tmp_box.expandBy(tmp_max);
    sub_boxes[1] = tmp_box;

    tmp_min = Point3f(half_pt.x(), min_pt.y(), half_pt.z());
    tmp_max = Point3f(max_pt.x(), half_pt.y(), max_pt.z());
    tmp_box.reset();
    tmp_box.expandBy(tmp_min);
    tmp_box.expandBy(tmp_max);
    sub_boxes[2] = tmp_box;

    tmp_min = Point3f(min_pt.x(), min_pt.y(), half_pt.z());
    tmp_max = Point3f(half_pt.x(), half_pt.y(), max_pt.z());
    tmp_box.reset();
    tmp_box.expandBy(tmp_min);
    tmp_box.expandBy(tmp_max);
    sub_boxes[3] = tmp_box;

    tmp_min = Point3f(min_pt.x(), half_pt.y(), min_pt.z());
    tmp_max = Point3f(half_pt.x(), max_pt.y(), half_pt.z());
    tmp_box.reset();
    tmp_box.expandBy(tmp_min);
    tmp_box.expandBy(tmp_max);
    sub_boxes[4] = tmp_box;

    tmp_min = Point3f(half_pt.x(), half_pt.y(), min_pt.z());
    tmp_max = Point3f(max_pt.x(), max_pt.y(), half_pt.z());
    tmp_box.reset();
    tmp_box.expandBy(tmp_min);
    tmp_box.expandBy(tmp_max);
    sub_boxes[5] = tmp_box;

    tmp_min = half_pt;
    tmp_max = max_pt;
    tmp_box.reset();
    tmp_box.expandBy(tmp_min);
    tmp_box.expandBy(tmp_max);
    sub_boxes[6] = tmp_box;

    tmp_min = Point3f(min_pt.x(), half_pt.y(), half_pt.z());
    tmp_max = Point3f(half_pt.x(), max_pt.y(), max_pt.z());
    tmp_box.reset();
    tmp_box.expandBy(tmp_min);
    tmp_box.expandBy(tmp_max);
    sub_boxes[7] = tmp_box;

    // 做8个TriArray
    TriArray tri_lists[8];
    int estimate_size = triangles.size() / 8;
    if (estimate_size > 0)
    {
        for (int i=0; i<8; i++)
        {
            tri_lists[i].reserve(estimate_size);
        }
    }
    // 判断是否相交
    for (auto tri_itr=triangles.begin(); tri_itr != triangles.end(); tri_itr++)
    {
        for (size_t i=0; i<8; i++)
        {
            BoundingBox3f cur_box = sub_boxes[i];
            if (this->overlaps(cur_box, *tri_itr, meshes))
            {
                tri_lists[i].push_back(*tri_itr);
            }
        }
    }

    // 生成子节点
    node_ptr->m_sub_nodes.reserve(8);
    for (int i=0; i<8; i++)
    {
        OctTreeNodePtr tmp_ptr = make_tree(sub_boxes[i], tri_lists[i], meshes);
        node_ptr->m_sub_nodes.push_back(tmp_ptr);
    }
    return node_ptr;
}

bool OctTree::overlaps(const BoundingBox3f &bounding_box, TriangleIndex triangle, const MeshArray &meshes)
{
    Mesh* mesh = meshes.at(triangle.mesh_idx);
    const MatrixXf &V  = mesh->getVertexPositions();
    const MatrixXu &F  = mesh->getIndices();
    uint32_t f = triangle.triangle_idx;
    uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);
    Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

    BoundingBox3f tri_box;
    tri_box.expandBy(p0);
    tri_box.expandBy(p1);
    tri_box.expandBy(p2);

    return bounding_box.overlaps(tri_box);
}

bool OctTree::nodeRayIntersect(OctTreeNodePtr node, Ray3f &ray, float &u, float &v, float &t, TriangleIndex &tri_idx) const
{
    bool foundIntersection = false;
    if (node == nullptr)
    {
        return foundIntersection;
    }


    if (!node->bounding_box.rayIntersect(ray))
        return false;
    
    if (node->is_leaf())
    {
        for (auto tri_itr=node->m_tri_arr.begin(); tri_itr != node->m_tri_arr.end(); tri_itr++)
        {
            int mesh_idx = tri_itr->mesh_idx;
            Mesh *mesh = m_mesh_array.at(mesh_idx);
            float tmp_u, tmp_v, tmp_t;
            bool is_intersect = mesh->rayIntersect(
                tri_itr->triangle_idx,
                ray,
                tmp_u,
                tmp_v,
                tmp_t
            );
            if (is_intersect)
            {
                u = tmp_u;
                v = tmp_v;
                t = tmp_t;
                ray.maxt = t;
                foundIntersection = true;
                tri_idx.mesh_idx = mesh_idx;
                tri_idx.triangle_idx = tri_itr->triangle_idx;
            }
        }
        return foundIntersection;
    }
    
    for (auto node_itr = node->m_sub_nodes.begin(); node_itr != node->m_sub_nodes.end(); node_itr++)
    {
        bool is_intersect = nodeRayIntersect(
            *node_itr,
            ray,
            u,
            v,
            t,
            tri_idx
        );
        if (is_intersect)
        {
            foundIntersection = is_intersect;
        }
    }
    return foundIntersection;
}

bool OctTree::rayIntersect(Ray3f &ray, float &u, float &v, float &t, TriangleIndex &tri_idx) const
{
    bool foundIntersection = false;
    foundIntersection = nodeRayIntersect(
        root,
        ray,
        u,
        v,
        t,
        tri_idx
    );
    fclose(m_log_file);
    return foundIntersection;
}

void Accel::addMesh(Mesh *mesh) {
    m_mesh = mesh;
    m_mesh_array.push_back(mesh);
    m_bbox.expandBy(mesh->getBoundingBox());
}

void Accel::build() {
    m_octtree.build(m_bbox, m_mesh_array);
}

// bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
//     bool foundIntersection = false;  // Was an intersection found so far?
//     uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

//     Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

//     /* Brute force search through all triangles */
//     for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
//         float u, v, t;
//         if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
//             /* An intersection was found! Can terminate
//                immediately if this is a shadow ray query */
//             if (shadowRay)
//                 return true;
//             ray.maxt = its.t = t;
//             its.uv = Point2f(u, v);
//             its.mesh = m_mesh;
//             f = idx;
//             foundIntersection = true;
//         }
//     }

//     if (foundIntersection) {
//         /* At this point, we now know that there is an intersection,
//            and we know the triangle index of the closest such intersection.

//            The following computes a number of additional properties which
//            characterize the intersection (normals, texture coordinates, etc..)
//         */

//         /* Find the barycentric coordinates */
//         Vector3f bary;
//         bary << 1-its.uv.sum(), its.uv;

//         /* References to all relevant mesh buffers */
//         const Mesh *mesh   = its.mesh;
//         const MatrixXf &V  = mesh->getVertexPositions();
//         const MatrixXf &N  = mesh->getVertexNormals();
//         const MatrixXf &UV = mesh->getVertexTexCoords();
//         const MatrixXu &F  = mesh->getIndices();

//         /* Vertex indices of the triangle */
//         uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

//         Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

//         /* Compute the intersection positon accurately
//            using barycentric coordinates */
//         its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

//         /* Compute proper texture coordinates if provided by the mesh */
//         if (UV.size() > 0)
//             its.uv = bary.x() * UV.col(idx0) +
//                 bary.y() * UV.col(idx1) +
//                 bary.z() * UV.col(idx2);

//         /* Compute the geometry frame */
//         its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

//         if (N.size() > 0) {
//             /* Compute the shading frame. Note that for simplicity,
//                the current implementation doesn't attempt to provide
//                tangents that are continuous across the surface. That
//                means that this code will need to be modified to be able
//                use anisotropic BRDFs, which need tangent continuity */

//             its.shFrame = Frame(
//                 (bary.x() * N.col(idx0) +
//                  bary.y() * N.col(idx1) +
//                  bary.z() * N.col(idx2)).normalized());
//         } else {
//             its.shFrame = its.geoFrame;
//         }
//     }

//     return foundIntersection;
// }

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    // traverse the octTree
//     for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
//         float u, v, t;
//         if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
//             /* An intersection was found! Can terminate
//                immediately if this is a shadow ray query */
//             if (shadowRay)
//                 return true;
//             ray.maxt = its.t = t;
//             its.uv = Point2f(u, v);
//             its.mesh = m_mesh;
//             f = idx;
//             foundIntersection = true;
//         }
//     }
    float u, v, t;
    TriangleIndex tri_idx;
    if (m_octtree.rayIntersect(ray, u, v, t, tri_idx))
    {
        its.t = ray.maxt;
        its.uv = Point2f(u, v);
        foundIntersection = true;
        its.mesh = m_mesh_array.at(tri_idx.mesh_idx);
        f = tri_idx.triangle_idx;
    }

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

