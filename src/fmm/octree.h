/*
* octree.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#ifndef OCTREE_H_
#define OCTREE_H_

#include <cassert>
#include <vector>
#include <deque>
#include <algorithm>
#include <map>
#include <complex>
#include <iostream>
#include "../common/math_vector.h"
#include "../common/charge.h"
#include "../common/octree_indexer.h"

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/iterator/filter_iterator.hpp>

// fwd decl.
template<typename NodeT, typename CType> class Octree;

template<class ContentType>
class Node {

public:

    typedef Node<ContentType> thisType;
    friend class Octree<thisType,ContentType>;
    typedef std::vector<ContentType*> ContentList;

    // Constructor for root node
    Node(const Vector& _centre, double _edge_length, const OctreeIndexer& idxer) : my_idx(idxer),
                                    my_level(my_idx.get_level()),
                                    centre(_centre),
                                    edge_length(_edge_length),
                                    has_children(false),
                                    is_leaf(true),
                                    is_shadow(false),
                                    is_deleted(false),
                                    cached_nb_size(-1)
    {
        parent_ptr = NULL;
        for (int ii=0; ii < 8; ++ii)
        {
            children[ii] = NULL;
        }
        descendents = 0;
        
        return;
    }

    virtual ~Node() {}

    // Constructor for child nodes, based on parent
    Node(Node& parent, OctreeIndexer idxer, unsigned short child_id) : my_idx(idxer),
                                                    my_level(my_idx.get_level()),
                                                    edge_length(parent.get_edge_length()/2.0),
                                                    has_children(false),
                                                    is_leaf(true),
                                                    is_shadow(false),
                                                    is_deleted(false),
                                                    contents(),
                                                    cached_nb_size(-1)
    {
        // stash pointer to parent
        parent_ptr = &parent;
        
        // child id runs from 1 to 8
        assert(child_id >= 1 && child_id <= 8);

        for (int ii=0; ii < 8; ++ii)
        {
            children[ii] = NULL;
        }

        // number of descendents
        descendents = 0;

        // figure out the centre
        const Vector& parent_centre = parent.get_centre();
        centre = parent_centre;

        unsigned short id=child_id;
        if (id <=4)
        {
            centre.z -= edge_length/2.0;
        }
        else
        {
            centre.z += edge_length/2.0;
            id -= 4;
        }

        if (id <= 2)
        {
            centre.y -= edge_length/2.0;
        }
        else
        {
            centre.y += edge_length/2.0;
            id -= 2;
        }

        if (id == 1)
        {
            centre.x -= edge_length/2.0;
        }
        else
        {
            assert(id == 2);
            centre.x += edge_length/2.0;
        }

        // now run through the contents of the parent,
        // add anything which needs to be added to this
        // octree node.
        for (typename ContentList::const_iterator it=parent.get_contents().begin(), end=parent.get_contents().end();
                it != end;
                ++it)
        {
            // get the child which should receive the item
            unsigned short intended_child_id = parent.get_child_for_xyz(**it);
            if (intended_child_id == child_id) { add(*it); }
        }

    }

    inline void delete_contents_recursive()
    {
        contents.clear();
        for (int child_id=1; child_id <= 8; ++child_id)
        {
            if (children[child_id-1] != NULL) {
                children[child_id-1]->delete_contents_recursive();
            }
        }
 
    }

    inline void append_children(std::vector<OctreeIndexer>& list)
    {
        for (int child_id=1; child_id <= 8; ++child_id)
        {
            if (children[child_id-1] != NULL) {
                children[child_id-1]->append_children(list);
            }
        }
        
        for (int child_id=1; child_id <= 8; ++child_id)
        {
            if (children[child_id-1] != NULL) {
                list.push_back(children[child_id-1]->get_idx());
            }

        }
        return;
    }
    
    inline unsigned short get_child_for_xyz(const Vector& pt) const
    {
        unsigned short child_id=1;
        if (pt.x > centre.x) {child_id += 1;}
        if (pt.y > centre.y) {child_id += 2;}
        if (pt.z > centre.z) {child_id += 4;}

        return child_id;
    }

    inline unsigned short get_id_within_parent() const
    {
        unsigned short id_within_parent = 1;
        if (my_idx.get_x_idx() % 2) {id_within_parent += 1;}
        if (my_idx.get_y_idx() % 2) {id_within_parent += 2;}
        if (my_idx.get_z_idx() % 2) {id_within_parent += 4;}
        return id_within_parent;
    }

    inline OctreeIndexer get_child_idx(unsigned short child_id) const
    {
        assert(child_id >=1 && child_id <= 8);
        return my_idx.get_child_idx(child_id);
    }

    inline const OctreeIndexer& get_idx() const { return my_idx; }
    //inline unsigned long& get_hashed_idx() const { return my_hash_idx; }

    inline unsigned short get_level() const { return my_level; }
    inline bool is_root() const { return get_level() == 0; }

    inline const Vector& get_centre() const { return centre; }
    inline double get_edge_length() const { return edge_length; }

    inline void add(ContentType* thing) {

        Vector delta = 2.0*(static_cast<const Vector&>(*thing) - centre) / edge_length;
        //std::cout << dynamic_cast<const Vector&>(*thing) << centre << " " << edge_length << " " << delta << "\n";

        // check the the thing belongs in this node
        //std::cout << Vector(*thing) << " " << centre << " " << delta
		//		  << " " << edge_length << std::endl;
        delta.x = fabs(delta.x);
        delta.y = fabs(delta.y);
        delta.z = fabs(delta.z);
        const double limit = 1.0 + 1e-9;
        assert(delta.x <= limit && delta.y <= limit && delta.z <= limit);

        contents.push_back(thing);

        return;

    }

    inline int nb_size() const { return cached_nb_size; }
    inline void set_cached_nb_size(int val) const { cached_nb_size = val; }
    inline int count_descendents() const { return descendents; }

    inline size_t size() const { return contents.size(); }
    inline const ContentList& get_contents() const { return contents; }
    inline ContentList& get_contents() { return contents; }
    //const ContentList& get_neighbourhood_contents() const { return neighbourhood_contents; }
    inline bool empty() const { return size() == 0; }

    // Has children is obvious enough- if the node has children
    // then this is true.
    inline bool hasChildren() const { return has_children; }
    inline void set_has_children() { has_children = true; }
    inline void unset_has_children() { has_children = false; }

    // a leaf node is not the same as a childless node!
    // A leaf is a node where the number of contents, and possibly the
    // number of contents in the entire neighbourhood, is within a
    // specified number.  However, because the tree is adaptive, in order
    // that all nodes have a complete neighbourhood, it is possible that
    // a neighbouring node requires this one to be subdivided, even though
    // from an evaluation point of view we would not traverse any further
    // down the tree.
    inline bool isLeaf() const { return is_leaf; }
    inline void unset_is_leaf() { is_leaf = false; }
    inline void set_is_leaf() { is_leaf = true; }

    // A shadow node is any leaf which is a child (or grandchild, or great-grandchild etc.)
    // of a leaf node.  (See comment above for definition of a leaf node).  Thus a shadow node
    // may or may not be childless; however it is definitely *not* a leaf.
    inline bool isShadow() const {  assert(is_shadow == false || is_leaf == false); return is_shadow; }
    inline void set_shadow() { is_shadow = true; }
    inline void unset_shadow() { is_shadow = false; }

    inline bool isDeleted() const { return is_deleted; }
    inline void set_deleted() { is_deleted = true; }
    inline void unset_deleted() { is_deleted = false; }

    inline bool check_pt_is_in_cube(const Vector& pt) const
    {
        double half_this_edge_length = edge_length / 2.0;
        return (pt.x >= (centre.x - half_this_edge_length)) && (pt.x <= (centre.x + half_this_edge_length)) &&
            (pt.y >= (centre.y - half_this_edge_length)) && (pt.y <= (centre.y + half_this_edge_length)) &&
            (pt.z >= (centre.z - half_this_edge_length)) && (pt.z <= (centre.z + half_this_edge_length));
    }

    // An octree node has a Vector position (centre of the cube)
    // It is indexed by an OctreeIndexer
    // Various data is hung from the node:
    // - a list of things in the node (or pointers to them)
    // - a list of data structures which hold, e.g. fast multipole method arrays

    template <typename some_int_type>
    inline Node& get_child(some_int_type child_id) const
    {
        
        // child_id runs from 1 to 8 (not 0 to 7)
        if (children[child_id-1] == NULL)
        {
            throw BadIndexer();
        }
        
        return *(children[child_id-1]);
    }
    
    mutable Node* children[8];
    Node* parent_ptr;
    int descendents;
    
protected:

    const OctreeIndexer my_idx;
    const unsigned short my_level;
    //const unsigned long my_hash_idx;
    Vector centre;
    const double edge_length;
    bool has_children;

    bool is_shadow;
    bool is_leaf;
    bool is_deleted;

    // the contents of this octree node
    ContentList contents;
    mutable int cached_nb_size;
    

};

template <class NodeType=Node<Vector>, class ContentType=Vector>
class Octree
{

public:

    typedef NodeType NodeT;
    typedef std::map<OctreeIndexer, NodeT*> NodeList;
    typedef std::map<unsigned short, NodeList*> LevelList;
    typedef typename NodeT::ContentList::iterator content_iterator;

    Octree() : built_neighbourhoods(false), strict_item_limit(true) {}

    Octree(unsigned int mean_items_in_nbhood, const Vector& centre, double edge_length, OctreeIndexer _root_idx=OctreeIndexer(0,0,0,0)) :
        root_idx(_root_idx),
        mean_items_in_neighbourhood(mean_items_in_nbhood),
        max_depth(OctreeIndexer::MAX_LEVELS),
        universe_edge_length(edge_length),
        built_neighbourhoods(false),
        strict_item_limit(true)
    {

        root_level = root_idx.get_level();

        // create zero'th level node list
        level_list[root_level] = new NodeList();

        // create the root node
        double root_edge_length = universe_edge_length;
        root_ptr = new NodeT(centre, root_edge_length, root_idx);

        // insert the root node
        NodeList& node_list = *(level_list[root_level]);
        node_list[root_ptr->get_idx()] = root_ptr;

    }

    virtual ~Octree() {

        const unsigned short top_level = get_top_level();
        const unsigned short bottom_level = get_bottom_level();

        // iterate over all nodes, and delete them
        for (unsigned short level=top_level; level <= bottom_level; ++level)
        {
            NodeList& nodes = get_node_list(level);
            for (typename NodeList::iterator it=nodes.begin(), end=nodes.end(); it != end; ++it)
            {
                NodeT* node_ptr = it->second;
                delete node_ptr;
            }
            nodes.clear();
            delete &(nodes);
        }
        level_list.clear();

        // the ref_counted_thing store contains boost::shared_ptr containing
        // heap-allocated copies of the things being stored.  Clear it and
        // they should all evaporate
        ref_counted_thing_store.clear();

    }

    // This function will forcible re-init the tree with a different neighbourhood size
    // Only does pointer-copying and re-allocating of new Octree nodes, so shouldn't 
    // alter the underlying objects.  Hopefully.
    void rejig(unsigned int new_mean_items_in_hood)
    {
        
        // reset the number of items in the neighbourhood
        if (mean_items_in_neighbourhood == new_mean_items_in_hood) {
            std::cerr << "Tried to rejig the tree but asked for same number of items in neighbourhood as before (=" << mean_items_in_neighbourhood << ") --> ignoring." << std::endl;
            return;
        }
        
        // set the new mean items per neighbourhood
        mean_items_in_neighbourhood = new_mean_items_in_hood;        
        
        const unsigned short top_level = get_top_level();
        const unsigned short bottom_level = get_bottom_level();

        // get pointer to root_node_list -- going to need it in a sec
        NodeList* root_node_list = &(get_node_list(top_level));
        
        // iterate over all nodes, and delete them - they're obsolete
        for (unsigned short level=top_level+1; level <= bottom_level; ++level)
        {
            NodeList& nodes = get_node_list(level);
            for (typename NodeList::iterator it=nodes.begin(), end=nodes.end(); it != end; ++it)
            {
                NodeT* node_ptr = it->second;
                delete node_ptr;
            }
            nodes.clear();
            delete &(nodes);
        }
        
        // need to clear expired level lists
        level_list.clear();
        
        // re-insert the root node list
        level_list[root_level] = root_node_list;
        
        // reset root node properties to be shiny and new
        root_ptr->unset_has_children();
        root_ptr->unset_shadow();
        root_ptr->set_is_leaf();
        root_ptr->set_cached_nb_size(0); // just in case start adding new stuff after rejig...
        
        typename NodeT::ContentList& everything = root_ptr->get_contents();
        if (everything.size() != 0 && ref_counted_thing_store.empty())
        {
            // if nothing in the ref-counted store then need to copy 
            // the list from root_ptr before re-inserting because the
            // user didn't use a reference-counting add function to
            // put stuff in the tree.
            typename NodeT::ContentList new_list;
            new_list.insert(new_list.begin(), everything.begin(), everything.end());
            
            everything.clear(); // nuke the contents of the root node

            // now re-add everything to the tree using the (dangerous) raw pointer add function
            // (dangerous because there's no reference counting going on at all)
            for (typename NodeT::ContentList::const_iterator it=new_list.begin(), end=new_list.end(); it != end; ++it)
            {
                add(*it);
            }
            
        }
        else
        {
            // double check that everything in the tree is also in the ref-counted store
            assert(ref_counted_thing_store.size() == everything.size());
            
            // ok can nuke the contents of the root node, and then add everything from the
            // ref_counted_thing_store, using raw pointers instead of re-copying the
            // shared pointers, otherwise we'll end up with duplicates in the 
            // ref_counted_thing_store.  And we wouldn't want that.
            everything.clear();
            
            for (typename std::vector<boost::shared_ptr<ContentType> >::const_iterator it=ref_counted_thing_store.begin(), end=ref_counted_thing_store.end(); it != end; ++it)
            {
                add(it->get()); // get the raw pointer from the boost::shared_ptr
            }
        }
        
        built_neighbourhoods = false;
        
        return;
       
    }

    inline void set_max_depth(unsigned short depth_limit)
    {
        if (depth_limit < OctreeIndexer::MAX_LEVELS) {
            max_depth = depth_limit;
        }
        return;
    }

    inline double get_edge_length() const { return universe_edge_length; }

    template <typename X>
    void add(const std::vector<X>& many)
    {
        // insert one by one as usual
        for (int ctr=0; ctr < many.size(); ++ctr)
        {
            add(many[ctr]);
        }
    }
    
#if 0
    // add a whole bunch of items to the tree
    void add(const std::vector<ContentType>& many)
    {
        // insert as normal - do it quick and dirty using raw
        // pointers (so no overhead due to copying)
        for (int ctr=0; ctr < many.size(); ++ctr)
        {
            add(const_cast<ContentType*>(&(many[ctr])));
        }

        // now loop over leaf nodes and make the entries in each 
        // node contiguous in memory -- this will be slow as we're 
        // doing a lot of copies.
        size_t running_ctr=0;
        std::vector< boost::shared_ptr<ContentType> > new_thing_store;
        new_thing_store.reserve(this->size());
        for (int level=get_bottom_level(); level >= get_top_level(); --level)
        {
            const NodeList& nodes = get_node_list(level);
            for (typename NodeList::const_iterator it=nodes.begin(), end=nodes.end(); it != end; ++it)
            {
                if (it->second->isLeaf())
                {
                    const typename NodeT::ContentList& contents = it->second->get_contents();
                    for (int ii=0; ii < contents.size(); ++ii)
                    {
                        new_thing_store.push_back( boost::shared_ptr<ContentType>(new ContentType(*(contents[ii])) )); // copy construct
                    }
                }
            }
        }

        // cascade delete through tree
        ref_counted_thing_store.clear();
        root_ptr->delete_contents_recursive();

        // reinsert - tree structure should exist so
        // this should be fairly fast.
        for (int ctr=0; ctr < new_thing_store.size(); ++ctr)
        {
            add(new_thing_store[ctr]);
        }

        return;

    }
#endif

    void add(const ContentType& thing)
    {
        // copy the thing
        boost::shared_ptr<ContentType> ref_counted_ptr(new ContentType(thing));
        ref_counted_thing_store.push_back( ref_counted_ptr );

        // add thing to the root node, which will recurse down
        add(ref_counted_ptr.get(), *root_ptr);

    }

    // DANGER: user must ensure the pointer passed in remains valid
    // for lifetime of the Octree. Otherwise segfaults this way lie...
    // Recommend using the add(boost::shared_ptr<ContentType> thing_shared_ptr)
    // function instead.
    void add(ContentType* thing_ptr)
    {
        // add thing to the root node, which will recurse down --
        // note that we have not made a safe copy of thing being stored
        // so better pray that the user keeps the thing_ptr in scope
        add(thing_ptr, *root_ptr);
    }

    void add(boost::shared_ptr<ContentType> thing_shared_ptr)
    {
        // if user passes in a shared_ptr then that's a good thing-
        // stash a copy of it in the ref_counted_thing_store and we'll
        // retain a reference count of the object for the lifetime of the
        // octree.
        ref_counted_thing_store.push_back(thing_shared_ptr);

        // add the raw pointer to the octree
        add(thing_shared_ptr.get(), *root_ptr);
    }

    inline size_t size() const
    {
        return root_ptr->size();
    }

    // An octree is a collection of Node objects, indexed by OctreeIndexer
    inline NodeT& get_node(const OctreeIndexer& idxer) const
    {
        const unsigned short level = idxer.get_level();
        NodeList& nodes_on_level = get_node_list(level);
        typename NodeList::iterator find_it = nodes_on_level.find(idxer);
        if (find_it == nodes_on_level.end()) { throw BadIndexer(); }
        return *(nodes_on_level.find(idxer)->second);
    }

    // Can also access by Vector -- returns the leaf node
    // containing that point
    inline const NodeT& get_node(const Vector& pt, int max_depth=-1) const
    {
        if (max_depth==-1) { max_depth = get_bottom_level(); }
        
        const NodeT* node = root_ptr;
        int depth=0;
        while (node->isLeaf()==false && depth < max_depth)
        {
            unsigned short child_id = node->get_child_for_xyz(pt);
            const NodeT* next_node_down;
            try {
                next_node_down = dynamic_cast<NodeT*>(&(node->get_child(child_id)));
                ++depth;
            }
            catch (BadIndexer) { break; }
            node = next_node_down;
        }
        return *node;
    }
    inline NodeT& get_node(const Vector& pt)
    {
        NodeT* node = root_ptr;
        while (node->isLeaf()==false)
        {
            unsigned short child_id = node->get_child_for_xyz(pt);
            NodeT* next_node_down;
            try {
                next_node_down = dynamic_cast<NodeT*>(&(node->get_child(child_id)));
            }
            catch (BadIndexer) { break; }
            node = next_node_down;
        }
        return *node;
    }

    inline boost::shared_ptr<typename NodeT::ContentList> get_neighbourhood_contents(const NodeT& node) const
    {

        boost::shared_ptr<typename NodeT::ContentList> neighbourhood_ptr(new typename NodeT::ContentList());
        assert(neighbourhood_ptr.get() != NULL);
        get_neighbourhood_contents(node, *neighbourhood_ptr);
        return neighbourhood_ptr;
    }

    inline void get_neighbourhood_contents(const NodeT& node, typename NodeT::ContentList& neighbourhood) const
    {
        if (node.nb_size() != -1) {
            neighbourhood.reserve(neighbourhood.size() + node.nb_size());
        }
        
        OctreeIndexer node_idxer = node.get_idx();
        for (unsigned short nb=0; nb < 27; ++nb)
        {
            try 
            {
                OctreeIndexer neighbour_idxer = node_idxer.get_neighbour_idxer(nb);
                //std::cout << "Getting idxer: " << idxer << std::endl;
                const NodeT& neighbour_node = get_node(neighbour_idxer);
                const typename NodeT::ContentList& contents = neighbour_node.get_contents();
                
                //std::cout << "node has contents: " << contents.size() << std::endl;
                neighbourhood.insert(neighbourhood.end(), contents.begin(), contents.end());
            }
            catch (BadIndexer) {}
        }
        
        if (node.nb_size() == -1) {
            node.set_cached_nb_size(neighbourhood.size());
        }

        return;
    }

//
//	void get_neighbourhood_contents(const NodeT& node, std::vector<ContentType*>& neighbourhood) const
//	{
//		for (unsigned short nb=0; nb < 27; ++nb)
//		{
//			try {
//				OctreeIndexer idxer = node.get_idx().get_neighbour_idxer(nb);
//				//std::cout << "Getting idxer: " << idxer << std::endl;
//				const typename NodeT::ContentList& contents = get_node(idxer).get_contents();
//				//std::cout << "node has contents: " << contents.size() << std::endl;
//				neighbourhood.reserve(neighbourhood.size() + contents.size());
//				neighbourhood.insert(neighbourhood.end(), contents.begin(), contents.end());
//			}
//			catch (BadIndexer) {}
//		}
//
//	}

#ifdef PREHYDROPHOBIC
    inline const ContentType& get_nearest(const Vector& where) const
#else
	static double defdist(const Vector& v, const ContentType& c)
			{ return (v - c).length2(); }
	// This should probably be templated, but that upsets backward compatibility
    inline const ContentType& get_nearest(
		const Vector& where,
		double (*distance)(const Vector&, const ContentType&) = defdist) const
#endif // PREHYDROPHOBIC
    {
        if (size() == 0) { throw std::exception(); }
        const NodeT* node = &(get_node(where));
        std::vector<ContentType*> neighbourhood;
        while (true)
        {
            get_neighbourhood_contents(*node, neighbourhood);
            if (neighbourhood.size() > 0) { break; }
            const NodeT& n = *node;
            OctreeIndexer parent_idx = node->get_idx().get_parent_idx();
            node = &(get_node(parent_idx));
        }

        // find nearest item in neighbourhood
        const ContentType* winner = *(neighbourhood.begin());
#ifdef PREHYDROPHOBIC
        double winning_dist = (where - *winner).length2();
#else
		double winning_dist = distance(where, *winner);
#endif // PREHYDROPHOBIC
        for (typename std::vector<ContentType*>::const_iterator it=neighbourhood.begin()+1, end=neighbourhood.end();
            it != end;
            ++it)
        {
#ifdef PREHYDROPHOBIC
            double dist = (where - **it).length2();
#else
			double dist = distance(where, **it);
#endif // PREHYDROPHOBIC
            if (dist < winning_dist)
            {
                winner = *it;
                winning_dist = dist;
            }
        }

        return *winner;
    }

    inline const NodeT& get_root_node() const { return *root_ptr; }
    inline NodeT& get_root_node() { return *root_ptr; }

    inline size_t get_num_nodes_on_level(unsigned short level) const
    {
        return get_node_list(level).size();
    }

    inline unsigned short get_top_level() const {
        return root_level;
    }
    inline unsigned short get_bottom_level() const {
        return get_top_level() + level_list.size() - 1;
    }

    inline size_t get_total_items() const 
    {
        return root_ptr->size();
    }

    inline size_t get_total_num_nodes() const {
        size_t total=0;
        for (unsigned short level=get_top_level(); level <= get_bottom_level(); ++level)
        {
            total += get_num_nodes_on_level(level);
        }
        return total;
    }

    NodeList& get_node_list(unsigned short level) const
    {
        typename LevelList::iterator find_it = level_list.find(level);
        assert(find_it != level_list.end());
        return *(find_it->second);
    }

    struct isNonEmptyChildless {
        inline bool operator()(const typename NodeList::value_type& pair) {
            return pair.second->hasChildren()==false && pair.second->empty()==false;
        }
    };

    typedef boost::filter_iterator<isNonEmptyChildless, typename NodeList::iterator> childless_iterator;
    typedef boost::filter_iterator<isNonEmptyChildless, typename NodeList::const_iterator> const_childless_iterator;

    inline content_iterator contents_begin() { return root_ptr->get_contents().begin(); }
    inline content_iterator contents_end() { return root_ptr->get_contents().end(); }

    inline childless_iterator childless_begin(unsigned short level)
    {
        NodeList& nodes = get_node_list(level);
        return boost::make_filter_iterator<isNonEmptyChildless>(nodes.begin(), nodes.end());
    }
    inline const_childless_iterator childless_begin(unsigned short level) const
    {
        const NodeList& nodes = get_node_list(level);
        return boost::make_filter_iterator<isNonEmptyChildless>(nodes.begin(), nodes.end());
    }
    inline childless_iterator childless_end(unsigned short level)
    {
        NodeList& nodes = get_node_list(level);
        return boost::make_filter_iterator<isNonEmptyChildless>(nodes.end(), nodes.end());
    }
    inline const_childless_iterator childless_end(unsigned short level) const
    {
        const NodeList& nodes = get_node_list(level);
        return boost::make_filter_iterator<isNonEmptyChildless>(nodes.end(), nodes.end());
    }

    void build_neighbourhoods()
    {
        assert(built_neighbourhoods == false);
        built_neighbourhoods = true;

        // iterate over all nodes and check that all nodes have complete neighbourhoods
        for (int level=get_top_level(); level <= get_bottom_level(); ++level)
        {
            std::vector<NodeT*> nodes_with_children;
            NodeList& nlst = get_node_list(level);
            for (typename NodeList::iterator node_it=nlst.begin(), node_end=nlst.end(); node_it != node_end; ++node_it)
            {
                NodeT* node = node_it->second;
                if (node->hasChildren() == true && node->isShadow() == false) {
                    nodes_with_children.push_back(node);
                }
            }

            for (typename std::vector<NodeT*>::iterator it=nodes_with_children.begin(), end=nodes_with_children.end(); it != end; ++it)
            {
                NodeT& node_with_children = **it;
                subdivide_neighbourhood(node_with_children);
            }

        }
    }

    // this will calculate the neighbourhood contents, ignoring the cached values
    // (so it should always be accurate -- use cached_count_neighbourhood_contents
    // for greater speed, but beware that the cached neighbourhood sizes might be
    // stale).
    inline size_t count_neighbourhood_contents(const NodeT& node) const
    {
        size_t retval = node.size();
        for (unsigned short nb=0; nb < 27; ++nb)
        {
            if (nb==13) { continue; } // already counted the node itself (see init. of retval above)
            
            try {
                OctreeIndexer idxer = node.get_idx().get_neighbour_idxer(nb);
                retval += get_node(idxer).size();
                
            }
            catch (BadIndexer) {}
        }
        
        return retval;
    }

    size_t calc_neighbourhood_interacts() const
    {
        size_t total_interacts=0;
        
        // the total product of <self_node_contents>*<neighbourhood_contents> for
        // all leaf nodes in the tree.  Corresponds to the number of explicit interactions
        // in the FMM near-field.
        for (int level=get_top_level(); level <= get_bottom_level() && level < max_depth; ++level)
        {
            for (typename NodeList::const_iterator node_it=get_node_list(level).begin(), node_end=get_node_list(level).end(); node_it != node_end; ++node_it)
            {
                const NodeT& node = (*node_it->second);
                if (node.isLeaf() && node.isShadow()==false && node.isDeleted()==false) {
                    long nb_size = cached_count_neighbourhood_contents(node);
                    total_interacts += nb_size*node.size();
                }
            }
        }
        
        return total_interacts;
    }
    
    struct NodeSorter
    {
        inline bool operator()(const NodeT* ptr1, const NodeT* ptr2) const
        {
            return ptr1->size() > ptr2->size();
        }
    };
    
    struct NodeNBSorter
    {
        inline bool operator()(const NodeT* ptr1, const NodeT* ptr2) const
        {
            return ptr1->nb_size() > ptr2->nb_size();
        }
    };
 
    void optimize_neighbourhoods()
    {
        // This function makes sure that the maximum amount of work per object in the tree 
        // (i.e. the explicit near-field interactions) is bounded according to 
        // mean_items_in_neighbourhood. Individual regions of the tree might slightly violate 
        // that constraint locally, to compensate for very sparse bits of the tree
        // where the constraint is massively under-met (i.e. the tree is over-divided).
        
        // Additional info:
        // In case you're curious about why this should be a good idea, it's because when
        // we subdivide a neighbourhood to reduce the explicit interactions, we *massively*
        // reduce the number of explicit interactions, so we undershoot our target amount
        // of explicit work- leading to an overly-zealously divided tree. If you plot the 
        // number of explicit interactions (or the number of interaction-list interactions
        // for that matter) as a function of number of objects in tree, you will observe
        // steps where the tree gets over-divided for a while, before filling up.

        // So anyway, this function tries to fix those 'steps' by allowing some local
        // neighbourhoods to not divide according to the strict rule: we iteratively
        // refine a fixed fraction of the 'fat' neighbourhoods until we hit the right
        // average for the whole tree.
        
        // must have already called build_neighbourhoods.
        assert(built_neighbourhoods == true);

        // total number of items we care about is the sum of all items in leaf nodes
        // not including shadows
        size_t total_items = 0;
        for (int level=get_top_level(); level <= get_bottom_level(); ++level)
        {
            for (typename NodeList::iterator node_it=get_node_list(level).begin(), node_end=get_node_list(level).end(); node_it != node_end; ++node_it)
            {
                NodeT* node = node_it->second;
                if (node->isShadow() == false && node->isLeaf() && node->isDeleted() == false) {
                    total_items += node->size();
                }
            }
        }
        
        // Note that the calc_neighbourhood_interacts function will iterate over all leaf
        // nodes and use the cached_calculate_neighbourhood_size function, which will set
        // a 'cached' neighbourhood sizes within each leaf node, so it doesn't need 
        // recalculating.  However it does also have the side-effect of leaving those cached
        // values in the nodes: so if you add stuff to the tree then re-call this function
        // it won't notice the new additions.  But then you should only be calling this
        // optimize function once anyway, so no biggy.
        double average_explicit_per_item = static_cast<double>(calc_neighbourhood_interacts()) / total_items;
        
        //std::cout << "Start while loop: average expl. per item: " << average_explicit_per_item << " (threshold: " << mean_items_in_neighbourhood << ")" << std::endl;
        while (floor(average_explicit_per_item) > mean_items_in_neighbourhood)
        {
            //std::cout << "Working: Average expl. per item: " << average_explicit_per_item << " (threshold: " << mean_items_in_neighbourhood << ")" << std::endl;
            
            // Get a list of all leaf nodes which need to be subdivided
            std::vector<NodeT*> sorted_nodes;
            for (int level=get_top_level(); level <= get_bottom_level(); ++level)
            {
                for (typename NodeList::iterator node_it=get_node_list(level).begin(), node_end=get_node_list(level).end(); node_it != node_end; ++node_it)
                {
                    NodeT* node = node_it->second;
                    if (node->isShadow() == false && node->isLeaf() && node->empty()==false) {
                        // only add to our sorted list if it's overweight-- i.e. too much in neighbourhood
                        if (cached_count_neighbourhood_contents(*node) > mean_items_in_neighbourhood)
                        {
                            sorted_nodes.push_back(node);
                        }
                    }
                }
            }
            
            // sort according to fatness - large neighbourhoods first
            // NB: this  sort puts them in descending order as our predicate is reversed from the usual
            std::sort(sorted_nodes.begin(), sorted_nodes.end(), NodeNBSorter());
            
            int fatty_ctr=0;
            int estimate_necessary_divisions = 1 + static_cast<int>(ceil(sorted_nodes.size() * static_cast<double>(average_explicit_per_item - mean_items_in_neighbourhood) / average_explicit_per_item));
            for (typename std::vector<NodeT*>::iterator leaf_it=sorted_nodes.begin(), leaf_end=sorted_nodes.end();
                 leaf_it != leaf_end; ++leaf_it)
            {
                NodeT& leaf = **leaf_it;
                leaf.unset_is_leaf();

                // if the fat leaf node already has children, they will be shadows, so need
                // to unset their shadowness, and make them into proper leaves
                if (leaf.hasChildren()) {
                    for (unsigned short ii=1; ii <= 8; ++ii)
                    {
                        try{
                            //OctreeIndexer child_idxer = leaf.get_child_idx(ii);
                            NodeT& child_node = dynamic_cast<NodeT&>(leaf.get_child(ii));
                            assert(child_node.isShadow());
                            assert(child_node.isLeaf() == false);
                            child_node.unset_shadow();
                            child_node.set_is_leaf();
                        }
                        catch (BadIndexer) {}
                    }
                }
                else
                {
                    subdivide(leaf);
                }

                // ensure that children of this node have a full neighbourhood too
                subdivide_neighbourhood(leaf);
                
                if ((++fatty_ctr) >= estimate_necessary_divisions) { break; }
            }
    
            // recount- then will test against threshold in while clause
            average_explicit_per_item = static_cast<double>(calc_neighbourhood_interacts()) / total_items;
            //std::cout << "looping - average expl. per item: " << average_explicit_per_item << " (threshold: " << mean_items_in_neighbourhood << ")" << std::endl;
            
        }
        //std::cout << "done" << std::endl;
        
        // achieved desired average neighbourhood size.  Job done.
        return;

    }

    inline unsigned int get_mean_items_in_neighbourhood() const { return mean_items_in_neighbourhood; }

    void linearize(unsigned short linearize_to)
    {
        // make octree linear down to the specified level
        for (int level=get_top_level(); level < linearize_to; ++level)
        {
            // geta list of the empties
            NodeList& node_list = get_node_list(level);
            for (typename NodeList::iterator node_it=node_list.begin(), node_end=node_list.end(); 
                 node_it != node_end; ++node_it)
            {
                NodeT* node = node_it->second;
                if (node->hasChildren() == false)
                {
                    subdivide(*node);
                    node->unset_is_leaf();
                }
            }
        }

    }

    inline int erase_children(NodeT& node) 
    {
#if 1        
        std::vector<OctreeIndexer> erase_list;
        node.append_children(erase_list);
        
        // now delete the list of nodes
        for (std::vector<OctreeIndexer>::const_iterator erase_it=erase_list.begin(), erase_end=erase_list.end(); erase_it != erase_end; ++erase_it)
        {
            NodeList& node_list = get_node_list(erase_it->get_level());
            delete node_list[*erase_it]; // don't forget to delete the heap-allocated memory!
            node_list.erase(*erase_it);
        }

        // clear child node pointers (now invalid)
        for (int ii=0; ii < 8; ++ii)
        {
            node.children[ii] = NULL;
        }
        node.unset_has_children();
        node.set_is_leaf();
        
        return erase_list.size();
#else   
        std::vector<OctreeIndexer> erase_list;
        node.append_children(erase_list);
        node.set_is_leaf();
        
        // now delete the list of nodes
        for (std::vector<OctreeIndexer>::const_iterator erase_it=erase_list.begin(), erase_end=erase_list.end(); erase_it != erase_end; ++erase_it)
        {
            NodeT& n = get_node(*erase_it);
            n.set_deleted();
            n.unset_is_leaf();
            n.set_shadow();
        }

        return erase_list.size();
#endif
    }
    
    void remove_empty_nodes()
    {
        // Remove empty nodes from the octree (optimization)
        for (int level=get_bottom_level(); level >= get_top_level(); --level)
        {
            // get a list of the empties
            NodeList& node_list = get_node_list(level);
            std::vector<OctreeIndexer> empties;
            for (typename NodeList::iterator node_it=node_list.begin(), node_end=node_list.end(); node_it != node_end; ++node_it)
            {
                NodeT* node = node_it->second;
				// Do not remove the root node as there is a pointer to it
                if (node->empty() == true && node != root_ptr) {
                    empties.push_back(node_it->first);
                    if (node->parent_ptr) {
                        node->parent_ptr->children[node->get_id_within_parent() - 1] = NULL;
                    }
                }
            }

            // now delete the empties from the map
            for (std::vector<OctreeIndexer>::iterator erase_it=empties.begin(), erase_end=empties.end(); erase_it != erase_end; ++erase_it)
            {
                delete node_list[*erase_it]; // don't forget to delete the heap-allocated memory!
                node_list.erase(*erase_it);
            }
        }

        return;

    }
    
    inline typename NodeT::ContentList& get_everything() const { return root_ptr->get_contents(); }
    inline void set_strict_divisions() { strict_item_limit = true; }
    inline void unset_strict_divisions() { strict_item_limit = false; }

protected:

    inline void add(ContentType* thing_ptr, NodeT& node)
    {
        //std::cout << "adding " << thing_ptr << " to " << node.get_idx() << " which already contains " << node.size() << " objects" << std::endl;

        // add thing to node
        node.add(thing_ptr);

        // check node is not overflowing
        if (node.hasChildren() == false)
        {
            if (strict_item_limit)
            {

                // check if necessary to subdivide
                // this is actually quite complicated because
                // the children might also need to subdivide
                // if it so happens that everything goes into the
                // same octant of the cube
                // Also our constraint is the size of the neighbourhood, not the 
                // absolute number of items in the node itself, so this is a coarse
                // strategy.
                std::deque<NodeT*> chk_queue;
                chk_queue.push_back(&node);

                while (!chk_queue.empty())
                {
                    NodeT* n = chk_queue.front();
                    chk_queue.pop_front();
                    if (n->size() > mean_items_in_neighbourhood && n->get_idx().get_level() < max_depth)
                    {
                        subdivide(*n);
                        n->unset_is_leaf();

                        // now add the children to the chk_queue
                        for (int child_id=1; child_id <=8; ++child_id)
                        {
                            NodeT& n_ref = dynamic_cast<NodeT&>(n->get_child(child_id));
                            chk_queue.push_back(&n_ref);
                        }
                    }
                }
            }
            else
            {
                // faster version
                // don't bother checking children
                // must use optimize_neighbourhoods after finished inserting items
                NodeT* n = &node;
                if (n->size() > mean_items_in_neighbourhood*4 && n->get_idx().get_level() < max_depth)
                {
                    subdivide(*n);
                    n->unset_is_leaf();
                }
            }
        }
        else
        {
            // pass on to child
            unsigned short child_id = node.get_child_for_xyz(*thing_ptr);
            NodeT& child = dynamic_cast<NodeT&>(node.get_child(child_id));
            add(thing_ptr, child);
        }

        return;
    }

    inline void subdivide(NodeT& node)
    {
        if (node.hasChildren()) { return; } // already subdivided
        node.set_has_children();

        unsigned short new_level = node.get_idx().get_level() + 1;
        assert(new_level > get_top_level());
        while (new_level > get_bottom_level())
        {
            assert(level_list.find(get_bottom_level()+1) == level_list.end());
            level_list[get_bottom_level()+1] = new NodeList();
        }
        NodeList& node_list = get_node_list(new_level);

        // create children
        for (unsigned short child_id=1; child_id <=8; ++child_id)
        {
            OctreeIndexer child_idxer = node.get_child_idx(child_id);
            NodeT* child_node = new NodeT(node, child_idxer, child_id);
            node.children[child_id-1] = child_node; // store pointer internally
            node_list[child_idxer] = child_node;
        }
        
        NodeT* n = &node;
        do
        {
            n->descendents += 8;
            n = dynamic_cast<NodeT*>(n->parent_ptr);
        } while (n);
        
        return;
    }

    void subdivide_neighbourhood(NodeT& node)
    {

        OctreeIndexer node_idxer = node.get_idx();
        for (unsigned short ii=0; ii < 27; ++ii)
        {
            try {
                OctreeIndexer idxer = node_idxer.get_neighbour_idxer(ii);
                NodeT& neighbour_node = get_node(idxer);
                if (neighbour_node.hasChildren()==false)
                {

                    // if the neighbour node has no children- then surely
                    // it is a leaf, or a shadow
                    assert(neighbour_node.isLeaf() || neighbour_node.isShadow());

                    subdivide(neighbour_node);

                    // label those new nodes as shadows- also they are *not* leaves!
                    // also: delete any empty shadows
                    for (unsigned short child_id=1; child_id <=8; ++child_id)
                    {
                        NodeT& new_child_node = dynamic_cast<NodeT&>(neighbour_node.get_child(child_id));
                        new_child_node.unset_is_leaf();
                        new_child_node.set_shadow();

//                         // if the child node is empty, delete it - don't need empty shadows
//                         if (new_child_node.empty()) {
//                             delete get_node_list(child_idxer.get_level())[child_idxer]; // delete heap allocated memory
//                             get_node_list(child_idxer.get_level()).erase(child_idxer);
//                         }

                    }
                }
                else
                {
                    for (unsigned short child_id=1; child_id <=8; ++child_id)
                    {
                        //OctreeIndexer child_idxer = idxer.get_child_idx(child_id);
                        NodeT& new_child_node = dynamic_cast<NodeT&>(neighbour_node.get_child(child_id));
                        new_child_node.unset_deleted();
                    }
                }
            }
            catch (BadIndexer) {}
        }
    }

    inline size_t cached_count_neighbourhood_contents(const NodeT& node) const
    {
        int retval = 0;
        retval = node.nb_size();
        
        // if the neighbourhood size has not been calculated for this node
        // then calculate it and then cache it.
        // WARNING: don't add anything to the tree after doing this or the
        // neighbourhood counts will be wrong.  So you should really only
        // use this function (and functions which call it) on a tree you are
        // done adding things to.
        if (retval == -1)
        {
            retval = static_cast<int>(count_neighbourhood_contents(node));
            node.set_cached_nb_size(retval);
        }
        
        return static_cast<size_t>(retval);
    }

    OctreeIndexer root_idx;
    unsigned short root_level;
    double universe_edge_length;
    unsigned int mean_items_in_neighbourhood;
    unsigned short max_depth;
    mutable LevelList level_list;
    NodeT* root_ptr;
    bool built_neighbourhoods;
    bool strict_item_limit;

    // a container holding shared_ptr to all heap allocated copies of objects -
    // will delete all the heap-allocated objects when octree is deleted
    std::vector< boost::shared_ptr<ContentType> > ref_counted_thing_store;

};

#endif /* OCTREE_H_ */
