package group.matsen.phylohmc.pll;

import com.sun.jna.Pointer;
import com.sun.jna.Structure;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

import java.util.Arrays;
import java.util.List;

public class pll_partition_t extends Structure {
    public int tips;
    public int clv_buffers;
    public int states;
    public int sites;
    public int rate_matrices;
    public int prob_matrices;
    public int rate_cats;
    public int scale_buffers;
    public int attributes;
    /** C type : const unsigned int* */
    public IntByReference map;
    /** vectorization options */
    public size_t alignment;
    public int states_padded;
    /** C type : double** */
    public PointerByReference clv;
    /** C type : double** */
    public PointerByReference pmatrix;
    /** C type : double* */
    public DoubleByReference rates;
    /** C type : double* */
    public DoubleByReference rate_weights;
    /** C type : double** */
    public PointerByReference subst_params;
    /** C type : unsigned int** */
    public PointerByReference scale_buffer;
    /** C type : double** */
    public PointerByReference frequencies;
    /** C type : double* */
    public DoubleByReference prop_invar;
    /** C type : int* */
    public IntByReference invariant;
    /** C type : unsigned int* */
    public IntByReference pattern_weights;
    /** C type : int* */
    public IntByReference eigen_decomp_valid;
    /** C type : double** */
    public PointerByReference eigenvecs;
    /** C type : double** */
    public PointerByReference inv_eigenvecs;
    /** C type : double** */
    public PointerByReference eigenvals;
    /** tip-tip precomputation data */
    public int maxstates;
    public int log2_maxstates;
    public int log2_states;
    public int log2_rates;
    /** C type : unsigned char** */
    public PointerByReference tipchars;
    /** C type : unsigned char* */
    public Pointer charmap;
    /** C type : double* */
    public DoubleByReference ttlookup;
    /** C type : unsigned int* */
    public IntByReference tipmap;

    {
        setAutoSynch(false);
    }

    @Override
    protected List getFieldOrder() {
        return Arrays.asList("tips", "clv_buffers", "states", "sites", "rate_matrices", "prob_matrices", "rate_cats", "scale_buffers", "attributes", "map", "alignment", "states_padded", "clv", "pmatrix", "rates", "rate_weights", "subst_params", "scale_buffer", "frequencies", "prop_invar", "invariant", "pattern_weights", "eigen_decomp_valid", "eigenvecs", "inv_eigenvecs", "eigenvals", "maxstates", "log2_maxstates", "log2_states", "log2_rates", "tipchars", "charmap", "ttlookup", "tipmap");
    }

}
