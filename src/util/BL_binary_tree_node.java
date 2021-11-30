package util;

public class BL_binary_tree_node
{

    BL_ranges_binary_tree element;
    int value;
    BL_binary_tree_node left = null;
    BL_binary_tree_node right = null;

    public BL_binary_tree_node(BL_ranges_binary_tree iu, int value)
    {
        this.element = iu;
        this.value = value;
    }
}
