import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.Stack;

public class Main {

    public static void main(String[] args) throws Exception{

        ArrayList<int[][]> matrixList = new ArrayList<>();

        try (Scanner scanner = new Scanner(new File("src/input"))){
            int cntInput = 0;
            while (scanner.hasNext()) {
                String firstLine = scanner.nextLine();
                String[] nm = firstLine.trim().split(" ");
                if (nm.length != 2) throw new IOException("input test case P length is invalid.");
                int n = Integer.parseInt(nm[0]);
                int m = Integer.parseInt(nm[1]);
                if (n == 0 && m == 0) break;
                int[][] matrix = new int[n][n];
                for (int i = 0; i < m; i++ ){
                    if (scanner.hasNext()){
                        String nl = scanner.nextLine();
                        String[] abc = nl.trim().split(" ");
                        if(abc.length != 3) throw new IOException("input (abc) length is invalid.");
                        int a = Integer.parseInt(abc[0]) - 1;
                        int b = Integer.parseInt(abc[1]) - 1;
                        int c = Integer.parseInt(abc[2]);
                        if (c == 1) {
                            matrix[a][b] = 1;
                        } else if (c == 2) {
                            matrix[a][b] = 1;
                            matrix[b][a] = 1;
                        } else {
                            throw new IOException("input c is valid");
                        }
                    } else {
                        throw new IOException("cannot found abc input.");
                    }
                }
                matrixList.add(matrix);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        // print adj matrix
        for (int[][] m : matrixList){
            printMatrix(m);
            System.out.println("-------------------------------");
        }

        // print and do the output
        for (int[][] m : matrixList){
            findAnswerAndAddLeastEdge(m);
        }

        System.out.println("\nAnswer Matrix");
        // print answer matrix
        for (int[][] m : matrixList){
            printMatrix(m);
            System.out.println("-------------------------------");
        }

    }

    public static void printMatrix(int[][] matrix){
        int n = matrix.length;
        if (matrix[0].length != n) throw new ArrayStoreException("matrix is invalid.");
        System.out.print("  ");
        for (int i = 0; i < n; i++){
            System.out.print(i + 1 + " ");
        }
        System.out.println();
        for (int i = 0; i < n; i++) {
            System.out.print(i + 1 + " ");
            for (int j = 0; j < n; j++){
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    public static void findAnswer(int[][] matrix){
        int n = matrix.length;
        if (matrix[0].length != n) throw new ArrayStoreException("matrix is invalid.");

        for(int i = 0; i < n; i++){

            Stack<Integer> route = new Stack<>();
            boolean[] isVisited = new boolean[n];
            route.push(i);
            addAllEdge(matrix, n, route, isVisited);

        }

        // check scc
        ArrayList<ArrayList<Integer>> sccList = new ArrayList<>();
        boolean[] isVisited = new boolean[n];
        for (int i = 0; i < n; i++) {
            if (!isVisited[i]) {
                ArrayList<Integer> scc = new ArrayList<>();
                if (matrix[i][i] == 1) {
                    for (int j = 0; j < n; j++) {

                        if (matrix[i][j] == 1 && matrix[j][i] == 1) {
                            isVisited[j] = true;
                            scc.add(j);
                        }
                    }
                } else {
                    scc.add(i);
                    isVisited[i] = true;
                }
                sccList.add(scc);
            }
        }

        System.out.println("numbers of scc : " + sccList.size());
        for (ArrayList<Integer> scc : sccList) {
            System.out.print("[ ");
            for (int i : scc) {
                System.out.print(i + 1 + " ");
            }
            System.out.println("]");
        }
        if (sccList.size() == 1){
            System.out.println("output : 1");
        } else {
         System.out.println("output : 0");
        }
        System.out.println();
    }

    public static void addAllEdge(int[][] matrix, int n, Stack<Integer> route, boolean[] isVisited){

        int i = route.peek();
        // find all next edge that connect ith edge
        for (int j = 0; j < n; j++) {
            if (matrix[i][j] == 1 && !isVisited[j]) {

                // add all vertex in route that all connect with jth edge
                for (int r : route){
                    matrix[r][j] = 1;
                }

                // DFS recursive
                isVisited[j] = true;
                route.push(j);
                addAllEdge(matrix, n, route, isVisited);
                route.pop();

            }
        }
    }

    public static void findAnswerAndAddLeastEdge(int[][] matrix) throws Exception {
        int n = matrix.length;
        if (matrix[0].length != n) throw new ArrayStoreException("matrix is invalid.");

        for(int i = 0; i < n; i++){

            Stack<Integer> route = new Stack<>();
            boolean[] isVisited = new boolean[n];
            route.push(i);
            addAllEdge(matrix, n, route, isVisited);

        }

        // check scc
        ArrayList<ArrayList<Integer>> sccList = new ArrayList<>();
        boolean[] isVisited = new boolean[n];
        for (int i = 0; i < n; i++) {
            if (!isVisited[i]) {
                ArrayList<Integer> scc = new ArrayList<>();
                if (matrix[i][i] == 1) {
                    for (int j = 0; j < n; j++) {

                        if (matrix[i][j] == 1 && matrix[j][i] == 1) {
                            isVisited[j] = true;
                            scc.add(j);
                        }
                    }
                } else {
                    scc.add(i);
                    isVisited[i] = true;
                }
                sccList.add(scc);
            }
        }

        System.out.println("numbers of scc : " + sccList.size());
        for (ArrayList<Integer> scc : sccList) {
            System.out.print("[ ");
            for (int i : scc) {
                System.out.print(i + 1 + " ");
            }
            System.out.println("]");
        }

        System.out.println("\nStarting add edge..");
        // add the least edge
        addLeastEdge(matrix, sccList, n);

        // update and finish
        findAnswer(matrix);
    }

    public static void addLeastEdge(int[][] matrix, ArrayList<ArrayList<Integer>> sccList, int n) throws Exception {

        // Add edge method:
        // 1. create the largest cycle
        // 2. scan and match non-in-degree and non-out-degree
        // 3. connect fragment of non-in-degree/non-out-degree to some vertex in the largest cycle

        // method 1 :

        // create stack for non-in-degree
        Stack<Integer> nonInDegreeStack = new Stack<>();
        Stack<Integer> nonOutDegreeStack = new Stack<>();
        // create candidate vertex for each scc
        int[] vertexs = new int[sccList.size()];

        // scan vertex that has not in-degree by use candidate vertex for each scc
        for (int s = 0; s < sccList.size(); s++) {
            ArrayList<Integer> scc = sccList.get(s);
            boolean hasInDegree = false;

            // use first vertex among scc
            int i = scc.get(0);
            // assign candidate vertex
            vertexs[s] = i;

            for (int j = 0; j < n; j++) {

                // skip vertex which in scc
                if (scc.contains(j)) {
                    continue;
                }
                // find in-degree
                if (matrix[j][i] == 1) {
                    hasInDegree = true;
                    break;
                }
            }

            // add scc that non-in-degree
            if (!hasInDegree) {
                nonInDegreeStack.push(i);
            }

        }

        // DEBUG
        //System.out.println(nonInDegreeStack.toString());


        // create path for make cycle

        // variable available
        // nonInDegreeStack - stack of non-in-degree
        // vertexs - candidate of each vertex in scc
        boolean[] path = new boolean[n]; // path shows what vertex is visited
        int lastVertex; // last vertex in path that don't have out-degree
        int targetVertex; // vertex that may be is connected by lastVertex
        int foundIsVisit; // find vertex that is visited if 0 = false, 1 = true
        int[] result; // use for result of findPath() [lastVertex, foundIsVisit] (is visited if 0 = false, 1 = true)

        // assign root
        int root = nonInDegreeStack.pop();
        path[root] = true;
        // find path for root
        result = findPath(matrix, vertexs, root, path);
        lastVertex = result[0];
        foundIsVisit = result[1];
        if (foundIsVisit == 1) throw new Exception("Impossible.");

        while (!nonInDegreeStack.isEmpty()) {
            targetVertex = nonInDegreeStack.pop();

            result = findPath(matrix, vertexs, targetVertex, path);

            // critical session - if found vertex that is visited, cannot use this path
            if (result[1] == 0) {   // not found vertex that is visited

                // connect lastVertex to targetVertex
                matrix[lastVertex][targetVertex] = 1;
                System.out.println("add : " + lastVertex + " to " + targetVertex);

                // assign new lastVertex
                lastVertex = result[0];
            } // if found, this non-in-degree has been pop and this path was marked by isVisited in findPath method
        }
        // stack of non-in-degree is empty, so we will connect to root for cycle happen
        matrix[lastVertex][root] = 1;
        System.out.println("add : " + lastVertex + " to " + root);
        // complete method 1 : create the largest cycle

        // method 2 :

        // update matrix
        for(int i = 0; i < n; i++){

            Stack<Integer> route = new Stack<>();
            boolean[] isVisited = new boolean[n];
            route.push(i);
            addAllEdge(matrix, n, route, isVisited);

        }

        // scan again to add non-in-degree and non-out-degree to stack
        for (ArrayList<Integer> scc : sccList) {
            boolean hasInDegree = false;
            boolean hasOutDegree = false;

            // use first vertex among scc
            int i = scc.get(0);

            for (int j = 0; j < n; j++) {

                // skip vertex which in scc
                if (scc.contains(j)) {
                    continue;
                }
                // find out-degree
                if (matrix[i][j] == 1) {
                    hasOutDegree = true;
                }
                // find in-degree
                if (matrix[j][i] == 1) {
                    hasInDegree = true;
                }
                if (hasInDegree && hasOutDegree) {
                    break;
                }
            }

            // add scc that non-in-degree and non-out-degree
            if (!hasInDegree) {
                nonInDegreeStack.push(i);
            }
            if (!hasOutDegree) {
                nonOutDegreeStack.push(i);
            }
        }

        // match non-in-degree and non-out-degree
        while (!nonInDegreeStack.isEmpty() && !nonOutDegreeStack.isEmpty()) {
            lastVertex = nonOutDegreeStack.pop();
            targetVertex = nonInDegreeStack.pop();

            matrix[lastVertex][targetVertex] = 1;
            System.out.println("add : " + lastVertex + " to " + targetVertex);
        }

        // method 3 :

        // fragment of non-in-degree
        while (!nonInDegreeStack.isEmpty()) {
            targetVertex = nonInDegreeStack.pop();

            matrix[root][targetVertex] = 1;
            System.out.println("add : " + root + " to " + targetVertex);
        }
        // fragment of non-out-degree
        while (!nonOutDegreeStack.isEmpty()) {
            lastVertex = nonOutDegreeStack.pop();

            matrix[lastVertex][root] = 1;
            System.out.println("add : " + lastVertex + " to " + root);
        }

    }

    // return [lastVertex, foundIsVisit] (is visited if 0 = false, 1 = true)
    public static int[] findPath(int[][] matrix, int[] vertexs, int i, boolean[] isVisited) {

        // vertex i / lastVertex
        int vi = i;
        // find vertex that is visited if 0 = false, 1 = true
        int foundIsVisit = 0;

        for (int vj : vertexs) {

            // vj : vertex j

            // skip same vertex
            if (vj == vi) {
                continue;
            }
            // find that vertex i connect with vertex j
            if (matrix[vi][vj] == 1) {
                if (!isVisited[vj]) {
                    isVisited[vj] = true;
                    // DFS
                    int[] result = findPath(matrix, vertexs, vj, isVisited);
                    vi = result[0];
                    foundIsVisit = result[1];
                } else {
                    // found visit = true
                    foundIsVisit = 1;
                }
                break;
            }
        }

        // cannot find path anymore
        return new int[]{vi,foundIsVisit};
    }
}
