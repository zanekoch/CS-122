/**
 *  GLOBAL ALIGNMENT PROBLEM
 *
 *  Find the highest-scoring alignment between two strings using a scoring matrix.
 *
 *  Given:  Two protein strings.
 *  Return: The maximum alignment score of these strings followed by an alignment achieving this maximum score. Use the BLOSUM62 scoring matrix and indel penalty sigma = 5. (If multiple alignments achieving the maximum score exist, you may return any one.)
 */

import java.util.*;
public class GlobalAlignmentProblem {
    public static String aminosAlpha = "ACDEFGHIKLMNPQRSTVWY";
    static int[][] blosum62 = {{4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2},{0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2},{-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3},{-1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2},{-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3},{0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3},{-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2},{-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1},{-1,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2},{-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1},{-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1},{-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2},{-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3},{-1,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1},{-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2},{1,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-2},{0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2},{0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1},{-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2},{-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7}};

    public static void main( String[] args ) {
        // DEFINE "s1" AND "s2"!!!
        String s1 = "";
        String s2 = "";
        // DEFINE "s1" AND "s2"!!!
        int indel = 5;
        String[] out = proteinGlobalAlignment(s1,s2,indel);
        System.out.println(out[0]);
        System.out.println(out[1]);
        System.out.println(out[2]);
    }

    public static String[] proteinGlobalAlignment( String s, String t, int indel ) {
        // Uses BLOSUM62. PASS indel AS POSITIVE INTEGER, NOT NEGATIVE!
        int[][] S = new int[s.length()+1][t.length()+1];
        int[][] opt = new int[s.length()+1][t.length()+1]; // 1 = right, 2 = down, 3 = diag
        S[0][0] = 0;
        for(int i = 1; i <= s.length(); ++i) {
            S[i][0] = S[i-1][0] - indel;
            opt[i][0] = 1;
        }
        for(int j = 1; j <= t.length(); ++j) {
            S[0][j] = S[0][j-1] - indel;
            opt[0][j] = 2;
        }
        for(int i = 1; i <= s.length(); ++i) {
            for(int j = 1; j <= t.length(); ++j) {
                int opt1 = S[i][j-1] - indel; // right
                int opt2 = S[i-1][j] - indel; // down
                int opt3 = S[i-1][j-1] + blosum62[aminosAlpha.indexOf(s.charAt(i-1))][aminosAlpha.indexOf(t.charAt(j-1))]; // diag
                S[i][j] = opt1;
                opt[i][j] = 1;
                if(opt2 > S[i][j]) {
                    S[i][j] = opt2;
                    opt[i][j] = 2;
                }
                if(opt3 > S[i][j]) {
                    S[i][j] = opt3;
                    opt[i][j] = 3;
                }
            }
        }
        String[] out = new String[3];
        out[0] = "" + S[s.length()][t.length()];
        out[1] = "";
        out[2] = "";
        int i = s.length();
        int j = t.length();
        while(i > 0 && j > 0) {
            if(opt[i][j] == 1) { // right
                out[1] = '-' + out[1];
                out[2] = t.charAt(j-- - 1) + out[2];
            }
            else if(opt[i][j] == 2) { // down
                out[1] = s.charAt(i-- - 1) + out[1];
                out[2] = '-' + out[2];
            }
            else if(opt[i][j] == 3) { // diag
                out[1] = s.charAt(i-- - 1) + out[1];
                out[2] = t.charAt(j-- - 1) + out[2];
            }
        }
        if(i > 0) {
            out[1] = s.substring(0,i) + out[1];
            String add = "";
            for(int x = 0; x < i; ++x) {
                add += '-';
            }
            out[2] = add + out[2];
        }
        if(j > 0) {
            out[2] = t.substring(0,j) + out[2];
            String add = "";
            for(int x = 0; x < j; ++x) {
                add += '-';
            }
            out[1] = add + out[1];
        }
        return out;
    }
}
