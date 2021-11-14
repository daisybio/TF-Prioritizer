package util;

import java.util.ArrayList;

public class ENSG_binary_tree {
    public ENSG_binary_tree_node root;

    public ENSG_binary_tree(ENSG_binary_tree_node root) {
        this.root = root;
    }

    public void add(int value, ENSG_ranges_binary_trees iu) {
        root = addRecursive(root, value, iu);
    }

    public ENSG_ranges_binary_trees containsNode(ENSG_ranges_binary_trees iu) {
        return containsNodeRecursive(root, iu, root);
    }

    private ENSG_binary_tree_node addRecursive(ENSG_binary_tree_node current, int value, ENSG_ranges_binary_trees iu) {
        if (current == null) {
            return new ENSG_binary_tree_node(iu, value);
        }

        if (value < current.value) {
            current.left = addRecursive(current.left, value, iu);
        } else if (value > current.value) {
            current.right = addRecursive(current.right, value, iu);
        } else {
            // value already exists
            return current;
        }

        return current;
    }

    private ENSG_ranges_binary_trees containsNodeRecursive(ENSG_binary_tree_node current, ENSG_ranges_binary_trees iu,
                                                           ENSG_binary_tree_node node_before) {
        if (current == null) {
            //return node_before.element;
            return null;
        }
        if (current.element.left_border <= iu.left_border && current.element.right_border >= iu.right_border &&
                current.element.left_border <= iu.right_border && current.element.right_border >= iu.left_border) {
            return current.element;
        }

        if (current.element.left_border >= iu.left_border && current.element.left_border >= iu.right_border &&
                current.element.right_border >= iu.right_border && current.element.right_border >= iu.left_border) {
            node_before = current;
            return containsNodeRecursive(current.left, iu, node_before);
        } else {
            node_before = current;
            return containsNodeRecursive(current.right, iu, node_before);
        }
    }
}
