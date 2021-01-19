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

#pragma once

#include <nori/mesh.h>
#include <vector>
#include <memory>

NORI_NAMESPACE_BEGIN

class OctTree;
struct TriangleIndex;
class OctTreeNode;

// mesh array 用于存储场景的meshes
using MeshArray = std::vector<Mesh *>;
using TriArray = std::vector<TriangleIndex>;
using OctTreeNodeArray = std::vector<std::shared_ptr<OctTreeNode>>;

/**
 * \brief TriangeIndex 分为 mesh_idx, triangle_idx
 * mesh_idx 索引对应的mesh。 triangle_idx 索引对应mesh的 triangle
 */
struct TriangleIndex {
    int mesh_idx;
    int  triangle_idx;

    TriangleIndex() {mesh_idx = -1; triangle_idx = -1;}
    TriangleIndex(int mesh_idx, int tri_idx)
    {
        this->mesh_idx = mesh_idx;
        this->triangle_idx = tri_idx;
    }
};

class OctTreeNode {
    // OctTreeNode当是中间节点时 m_tri_arr是空的，m_sub_nodes非空
    // OctTreeNode是叶子节点时： m_tri_arr非空， m_sub_nodes为空，没有子节点
    // boundingbox是每个节点都需要的用于判断是否和ray相交所以每个节点都会有
    TriArray m_tri_arr;
    OctTreeNodeArray m_sub_nodes;
    BoundingBox3f bounding_box;
public:
    bool is_leaf(){return m_sub_nodes.size() == 0;}
    friend OctTree;
};

using OctTreeNodePtr = std::shared_ptr<OctTreeNode>;

class OctTree {
    OctTreeNodePtr root;
    MeshArray m_mesh_array;
public:
    void build(const BoundingBox3f& bounding_box, const MeshArray &meshes);
    OctTreeNodePtr make_tree(const BoundingBox3f& bounding_box, TriArray triangles, const MeshArray &meshes);
    // 返回是否ray和对应的triangle相交
    bool overlaps(const BoundingBox3f &bounding_box, TriangleIndex triangle, const MeshArray &meshes);
    bool rayIntersect(Ray3f &ray, float &u, float &v, float &t, TriangleIndex &tri_idx) const;
    bool nodeRayIntersect(OctTreeNodePtr node, Ray3f &ray, float &u, float &v, float &t, TriangleIndex &tri_idx) const;
};

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {
public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure (currently a no-op)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

private:
    MeshArray     m_mesh_array;
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
    OctTree m_octtree;
};

NORI_NAMESPACE_END
