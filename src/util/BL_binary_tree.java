package util;

public class BL_binary_tree {

    public BL_binary_tree_node root;

    public BL_binary_tree(BL_binary_tree_node root)
    {
        this.root = root;
    }

    public void add(int value, BL_ranges_binary_tree iu)
    {
        root = addRecursive(root, value, iu);
    }

    public BL_ranges_binary_tree containsNode( BL_ranges_binary_tree iu) {
        return containsNodeRecursive(root, iu);
    }

    private BL_binary_tree_node addRecursive(BL_binary_tree_node current, int value,BL_ranges_binary_tree iu)
    {
        if (current == null) {
            return new BL_binary_tree_node(iu,value);
        }

        if (value < current.value) {
            current.left = addRecursive(current.left, value,iu);
        } else if (value > current.value) {
            current.right = addRecursive(current.right, value,iu);
        } else {
            // value already exists
            return current;
        }

        return current;
    }

    private BL_ranges_binary_tree containsNodeRecursive(BL_binary_tree_node current,BL_ranges_binary_tree iu)
    {
        if (current == null) {
            return null;
        }

        boolean is_inside = false;
        if(current.element.left_border <= iu.left_border && current.element.right_border >= iu.right_border && current.element.left_border <= iu.right_border && current.element.right_border >= iu.left_border)
        {
            is_inside=true;
        }
        boolean overlap_right = false;
        if(current.element.left_border <= iu.left_border && current.element.right_border <= iu.right_border && current.element.left_border<=iu.right_border && current.element.right_border>=iu.left_border)
        {
            overlap_right=true;
        }
        boolean overlap_left = false;
        if(current.element.left_border >= iu.left_border && current.element.right_border >= iu.right_border && current.element.left_border<=iu.right_border && current.element.right_border>=iu.left_border)
        {
            overlap_left=true;
        }
        boolean bl_inside = false;
        if(current.element.left_border>= iu.left_border && current.element.right_border <= iu.right_border && current.element.left_border <= iu.right_border && current.element.right_border >= iu.left_border)
        {
            bl_inside=true;
        }

        if(is_inside || overlap_left || overlap_right || bl_inside)
        {
            return current.element;
        }

        if(current.element.left_border >= iu.left_border)
        {
            return containsNodeRecursive(current.left,iu);
        }
        else
        {
            return containsNodeRecursive(current.right,iu);
        }
    }


}
