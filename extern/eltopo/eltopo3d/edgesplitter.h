// ---------------------------------------------------------
//
//  edgesplitter.h
//  Tyson Brochu 2011
//  
//  Functions supporting the "edge split" operation: subdividing an edge into two shorter edges.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_EDGESPLITTER_H
#define EL_TOPO_EDGESPLITTER_H

// ---------------------------------------------------------
//  Nested includes
// ---------------------------------------------------------

#include <cstddef>
#include <vector>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class SurfTrack;
template<unsigned int N, class T> struct Vec;
typedef Vec<3,double> Vec3d;
typedef Vec<2,size_t> Vec2st;
typedef Vec<3,size_t> Vec3st;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Edge splitter object.  Splits "long" edges by introducing a new vertex at the midpoint, optionally offsetting for curvature 
/// preservation.
///
// ---------------------------------------------------------

class EdgeSplitter
{
    
public:
    
    /// Constructor
    ///
    EdgeSplitter( SurfTrack& surf, bool use_curvature, double max_curvature_multiplier );
    
    /// Split all long edges
    ///
    bool split_pass();
    
    /// Split edges opposite large angles
    ///
    bool large_angle_split_pass();
    

    /// Maximum edge length.  Edges longer than this will be subdivided.
    ///
    double m_max_edge_length;   

    /// Whether to scale by curvature when computing edge lengths, in order to refine high-curvature regions
    ///
    bool m_use_curvature;
    
    /// The maximum curvature scaling allowed
    ///
    double m_max_curvature_multiplier;
        
private:
    
    /// The mesh this object operates on
    /// 
    SurfTrack& m_surf;   
    
    /// Check collisions between the edge [neighbour, new] and the given edge 
    ///
    bool split_edge_edge_collision(size_t neighbour_index, 
                                   const Vec3d& new_vertex_position, 
                                   const Vec3d& new_vertex_smooth_position, 
                                   const Vec2st& edge );
    
    /// Determine if the new vertex introduced by the edge split has a collision along its pseudo-trajectory.
    ///
    bool split_triangle_vertex_collision( const Vec3st& triangle_indices, 
                                         const Vec3d& new_vertex_position, 
                                         const Vec3d& new_vertex_smooth_position, 
                                         size_t overlapping_vert_index, 
                                         const Vec3d& vert );
    
    /// Determine if the pseudo-trajectory of the new vertex has a collision with the existing mesh.
    ///
    bool split_edge_pseudo_motion_introduces_intersection( const Vec3d& new_vertex_position, 
                                                          const Vec3d& new_vertex_smooth_position, 
                                                          size_t edge,
                                                          size_t tri0, 
                                                          size_t tri1,
                                                          size_t vertex_a,
                                                          size_t vertex_b,
                                                          size_t vertex_c,
                                                          size_t vertex_d );
    /// Determine if edge should be allowed to be split
    ///    
    bool edge_is_splittable( size_t edge_index );
    
    /// Split an edge, using subdivision_scheme to determine the new vertex location, if safe to do so.
    ///
    bool split_edge( size_t edge );
    
};


#endif
