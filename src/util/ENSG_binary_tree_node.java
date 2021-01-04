package util;

public class ENSG_binary_tree_node {
    ENSG_ranges_binary_trees element;
    int value;
    ENSG_binary_tree_node left=null;
    ENSG_binary_tree_node right=null;

    public ENSG_binary_tree_node(ENSG_ranges_binary_trees iu, int value)
    {
        this.element=iu;
        this.value=value;
    }
}
