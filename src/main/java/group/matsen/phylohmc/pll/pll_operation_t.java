package group.matsen.phylohmc.pll;

import com.sun.jna.Structure;

import java.util.Arrays;
import java.util.List;

public class pll_operation_t extends Structure implements Structure.ByReference {

    public int parent_clv_index = 0;
    public int parent_scaler_index = -1;
    public int child1_clv_index = 0;
    public int child1_matrix_index = 0;
    public int child1_scaler_index = -1;
    public int child2_clv_index = 0;
    public int child2_matrix_index = 0;
    public int child2_scaler_index = -1;

    {
        setAutoSynch(false);
        setAutoWrite(true);
    }

    @Override
    protected List getFieldOrder() {
        return Arrays.asList("parent_clv_index", "parent_scaler_index", "child1_clv_index", "child1_matrix_index", "child1_scaler_index", "child2_clv_index", "child2_matrix_index", "child2_scaler_index");
    }
}
