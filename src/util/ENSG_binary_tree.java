package util;

public class ENSG_binary_tree
{
    public ENSG_binary_tree_node root;

    public ENSG_binary_tree(ENSG_binary_tree_node root)
    {
        this.root = root;
    }

    public void add(int value, ENSG_ranges_binary_trees iu) {
        root = addRecursive(root, value, iu);
    }

    private ENSG_binary_tree_node addRecursive(ENSG_binary_tree_node current, int value,ENSG_ranges_binary_trees iu)
    {
        if (current == null) {
            return new ENSG_binary_tree_node(iu,value);
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
}
