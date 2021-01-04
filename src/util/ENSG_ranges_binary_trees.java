package util;

import weka.classifiers.trees.j48.EntropyBasedSplitCrit;

import java.util.HashSet;

public class ENSG_ranges_binary_trees implements Comparable
{
    public String chromosome ="";
    public int number = 0;
    public int left_border = 0;
    public int right_border = 0;
    public HashSet<String> ensgs = new HashSet<>();

    public boolean lies_between(Object o1, Object o2)
    {
        boolean isSame=false;
        ENSG_ranges_binary_trees other1 = ((ENSG_ranges_binary_trees)o1);
        ENSG_ranges_binary_trees other2 = ((ENSG_ranges_binary_trees)o2);

        if(other1.right_border <= this.left_border && other2.left_border >= this.left_border)
        {
            isSame=true;
        }

        return isSame;
    }

    public boolean is_match(Object o)
    {
        boolean isSame=false;
        ENSG_ranges_binary_trees other = ((ENSG_ranges_binary_trees)o);

        if(other.left_border==this.left_border&&other.right_border==this.right_border)
        {
            isSame=true;
        }

        if(this.left_border <= other.right_border && this.right_border >= other.left_border && this.left_border <= other.left_border && this.right_border <= other.right_border)
        {
            isSame=true;
        }

        if(this.left_border <= other.right_border && this.right_border >= other.left_border && this.left_border >= other.left_border && this.right_border >= other.right_border)
        {
            isSame=true;
        }

        if(this.left_border <= other.right_border && this.right_border >= other.left_border && this.left_border <= other.left_border && this.right_border>= other.right_border)
        {
            isSame=true;
        }



        return isSame;
    }

    public boolean findENSG(String ensg)
    {
        if(ensgs.contains(ensg))
        {
            return true;
        }
        return false;
    }

    public boolean isSameRange(Object o)
    {
        boolean isSame=false;
        ENSG_ranges_binary_trees other = ((ENSG_ranges_binary_trees)o);

        if(other.left_border==this.left_border && other.right_border==this.right_border)
        {
            isSame=true;
        }
        return isSame;
    }

    public  boolean isOverLap(Object o)
    {
        boolean isSame=false;
        ENSG_ranges_binary_trees other = ((ENSG_ranges_binary_trees)o);

        if(other.left_border < this.right_border)
        {
            isSame=true;
        }


        return isSame;
    }

    public String toString_binary()
    {
        String ret = "";

        ret += number + "\t" + left_border + "\t" + right_border + "\t";

        int count = 0;
        for(String s: ensgs)
        {
            if(count>0)
            {
                ret += ";" + s;
            }
            else
            {
                ret+=s;
            }

            count++;
        }

        return ret;
    }
    public String toString()
    {
        String ret = "";

        ret += chromosome + "\t" + left_border + "\t" + right_border + "\t";

        int count = 0;
        for(String s: ensgs)
        {
            if(count>0)
            {
                ret += ";" + s;
            }
            else
            {
                ret+=s;
            }

            count++;
        }

        return ret;
    }

    public boolean isTheSame(Object o)
    {
        boolean isSame = false;

        ENSG_ranges_binary_trees other = ((ENSG_ranges_binary_trees)o);

        if(this.left_border==other.left_border && this.right_border== other.right_border)
        {
            boolean same_ensgs_1 = true;
            boolean same_ensgs_2 = true;

            for(String s: this.ensgs)
            {
                if(!other.ensgs.contains(s))
                {
                    same_ensgs_1=false;
                }
            }

            for(String s: other.ensgs)
            {
                if(!this.ensgs.contains(s))
                {
                    same_ensgs_2=false;
                }
            }

            if(same_ensgs_1&&same_ensgs_2)
            {
                isSame=true;
            }

        }


        return isSame;

    }

    public String ensgs_to_String()
    {
        String ret = "";

        int count =0;
        for(String s: ensgs)
        {

            if(count>0)
            {
                ret+=";"+s;
            }
            else
            {
                ret+=s;
            }

            count++;
        }

        return ret;
    }

    @Override
    public int compareTo(Object o) {
        int compare = ((ENSG_ranges_binary_trees)o).left_border;

        return this.left_border-compare;
    }
}
