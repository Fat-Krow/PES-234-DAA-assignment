#include "header.h"
#include <limits.h>
static void reverse_graph(int n, const connection_t connection[n][n], connection_t rev[n][n]); // returns a transposed digraph i.e reverses the edges in the digraph
static void dfs(int i, int n, const connection_t connection[n][n], int airports[n], int stack[n], int *top);
static int count_ones(int n, int arr[n]);
static int path_find(int n, int dest, const connection_t connection[n][n], int k, int sour, int *visited, int curr_path_length);
static void swap(airport_t *a, airport_t *b);
static int partition(int n, airport_t array[n], int (*predicate_func)(const airport_t *, const airport_t *), int low, int high);
static void quickSort(int n, airport_t array[n], int (*predicate_func)(const airport_t *, const airport_t *), int low, int high);
static void definently_not_strncpy(char *dest, const char *src, int k);
static int definently_not_strlen(const char *string);
static int binarySearchCount(int n, const int arr[n], int key);
static void shifttable(int shift_table[256], const char pattern[]);
static int horspool(const char src[], const char p[]);
typedef struct
{
    int cost;
    int sequence_of_visits[8];
} trip;
static int get_all_trips(int n,
                         const connection_t connections[n][n],
                         int visited[],
                         int start,
                         int currPos,
                         int exempted,
                         int visit_count,
                         int cost,
                         trip all_trips[n],
                         int *size);
static void extract_trip_seq(int n,
                             int visited[n],   // input array
                             int trip_seq[n]); // output array
static int min(int n, trip all_trips[n]);

typedef struct edge_maker
{
    int u;
    int v;
    connection_t w;
} edge;

typedef struct edge_list_maker
{
    edge data[100];
    int n;
} edge_list;
static void bubble_sort(edge_list *The_list_of_all_edges);
static int find(int parent_track[], int vertexno);
static void applyUnion(int n, int parent_track[], int c1, int c2);
static void Kruskals_Algorithim(int n, pair_t edges[n - 1], const connection_t connections[n][n], edge_list *The_list_of_all_edges, edge_list *The_span_list);
static int min(int n, trip all_trips[n]);
static int min_two(int x, int y);
// ANY STATIC FUNCTIONS ARE UP HERE
static void reverse_graph(int n, const connection_t connection[n][n], connection_t rev[n][n])
{ // connection_t copy;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            rev[i][j] = connection[j][i];
        }
    }
}
static void dfs(int i, int n, const connection_t connection[n][n], int airports[n], int stack[n], int *top)
{
    airports[i] = 1;

    for (int j = 0; j < n; j++)
    { // edge=connection[i][j];
        if ((connection[i][j].distance != INT_MAX) && (connection[i][j].time != INT_MAX))
        {
            if (airports[j] == 0)
            {
                dfs(j, n, connection, airports, stack, top);
            }
        }
    }

    stack[*top] = i;
    *top = *top + 1;
}
int count_ones(int n, int arr[n])
{
    int count = 0;
    for (int i = 0; i < n; i++)
    {
        if (arr[i] == 1)
        {
            count = count + 1;
        }
    }
    return count;
}
static int path_find(int n, int dest, const connection_t connection[n][n], int k, int src, int *visited, int curr_path_length)
{
    if (src == dest)
    {
        return 1; // if src== dest we have found the way
    }

    if (k < curr_path_length)
    {
        return 0; // test path faild
    }
    visited[src] = 1;

    for (int j = 0; j < n; j++)
    {
        if (visited[j] != 1 && connection[src][j].distance != INT_MAX && connection[src][j].time != INT_MAX)
        {
            int flag = path_find(n, dest, connection, k, j, visited, curr_path_length + 1);
            // recursively call path find to see if a path exists from the adjacent nodes
            if (flag == 1)
            {
                return 1; // we found the path
            }
        }
    }
    return 0; // from this node there is no unvisted node
}
// function to swap elements
static void swap(airport_t *a, airport_t *b)
{
    airport_t t = *a;
    *a = *b;
    *b = t;
}
// function to find the partition position
static int partition(int n, airport_t array[n], int (*predicate_func)(const airport_t *, const airport_t *), int low, int high)
{

    // first I select the rightmost element as pivot
    airport_t pivot = array[high];

    // this is a  pointer for greater element
    int i = (low - 1);

    // then I  traverse each element of the array
    //&& compare them with the pivot
    for (int j = low; j < high; j++)
    {
        if (predicate_func(&array[j], &pivot))
        {

            // if a element < pivot is found
            // swap it with the >er element pointed by i
            i++;

            // swaping arr[i]  with arr[j]
            swap(&array[i], &array[j]);
        }
    }

    // and then finally I  swap the pivot element with the greater element located at  i
    swap(&array[i + 1], &array[high]);

    //&& return the partition point
    return (i + 1);
}
static void quickSort(int n, airport_t array[n], int (*predicate_func)(const airport_t *, const airport_t *), int low, int high)
{

    if ((low < high) && predicate_func(&array[low], &array[high])) // this is for ascending
    {

        // first I  find the pivot element such that
        //  elements < pivot go  left of pivot
        //&& elements > than pivot are on right
        int pi = partition(n, array, predicate_func, low, high);

        // left side call

        quickSort(n, array, predicate_func, low, pi - 1);

        // right side call
        quickSort(n, array, predicate_func, pi + 1, high);
    }
    else if ((low < high) && predicate_func(&array[high], &array[low])) // this for descending patterns
    {                                                                   // this low<high is required to avoid a outofbound error

        // first I  find the pivot element such that
        //  elements < pivot go  left of pivot
        //&& elements > than pivot are on right
        int pi = partition(n, array, predicate_func, low, high);

        // left side call

        quickSort(n, array, predicate_func, low, pi - 1);

        // right side call
        quickSort(n, array, predicate_func, pi + 1, high);
    }
}
static void definently_not_strncpy(char *dest, const char *src, int k)
{
    int i;
    for (i = 0; ((i < k) && (src[i] != '\0')); i++)
    {
        dest[i] = src[i]; // definently_not_strncpy indeed
    }
    dest[i] = '\0';
}

// Expects the input array to contain consecutive integers in the range
// [1, n], but scattered in some arbitrary order.
// The successive elements of the output array when applied as indexes
// upon the input array, will give a sorted version of the input array.
static void extract_trip_seq(int n,
                             int visited[n],  // input array
                             int trip_seq[n]) // output array
{

    for (int i = 0; i <= n; i++)
    {
        // Look for the value i, in the input array.
        int j;
        int must_continue;
        for (j = 0, must_continue = 1; (j < n) && must_continue; )
        {
            if (visited[j] == i)
            {
                must_continue = 0;
            }
            else
            {
                j += 1;
            }
        }

        trip_seq[i] = j;
    }
} // end of extract_trip_seq()

static int definently_not_strlen(const char *string)
{
    unsigned int count = 0;
    while (*string != '\0')
    {
        count++;
        string++; ////definently_not_strlen indeed
    }
    return count;
}
// A modifed  binary search approach  to return
// the number of elements <=given key
static int binarySearchCount(int n, const int arr[n], int key)
{

    int left = 0;
    int right = n - 1;

    int count = 0; // this var will store the number of elements <= the key
    int mid;
    while (left <= right)
    {
        mid = (right + left) / 2;

        // Check if middle element is
        // less than or equal to key
        if (arr[mid] <= key)
        {
            /*basicaly there should be At least (mid + 1) elements are there
             that  are <=key*/
            count = mid + 1;
            left = mid + 1;
        }

        // If key is small, ignore right side
        else
        {
            right = mid - 1;
        }
    }

    return count;
}
static void shifttable(int shift_table[256], const char pattern[])
{ /*this function will create the shift table */
    int i, j, m;
    m = definently_not_strlen(pattern); // m=length of pattern
    for (i = 0; i < 256; i++)
    {
        shift_table[i] = m; // by default filling up the table with the length of pattern string
    }
    for (j = 0; j < m - 1; j++)
    {
        shift_table[(int)pattern[j]] = m - 1 - j;
    }
}
static int horspool(const char text[], const char pattern[])
{
    int shift_table[256];
    shifttable(shift_table, pattern);
    int i, k, m, n;
    n = definently_not_strlen(text);    // length of the string IN WHICH we need to search
    m = definently_not_strlen(pattern); // length of the pattern we need to search for
    i = m - 1;                          // because we start from last char
    // int where;
    // shifttable(shift_table,pattern);//this function prepares the shift table
    while (i < n)
    {
        k = 0; // first char of pattern string

        while ((k < m) && (pattern[m - 1 - k] == text[i - k])) // i=2-1
        {
            k = k + 1;
            /*the real horspool logic starts here by comparing the END to BEGINN of PATTERN string
            with corresponding text string char
            */
        }
        if (k == m)
        {

            return (i - m + 1); // returns where it is found
        }
        else
        { /*This one line is makes horspool diffrent from any brute force
          comparisson tequique*/
            i = i + shift_table[(int)text[i]];
        }
    }
    return -1;
}

static void bubble_sort(edge_list *The_list_of_all_edges)
{
    int i, j;
    edge temp;

    for (i = 1; i < The_list_of_all_edges->n; i++)
    {
        for (j = 0; j < The_list_of_all_edges->n - 1; j++)
        {
            if (The_list_of_all_edges->data[j].w.time > The_list_of_all_edges->data[j + 1].w.time)
            {
                temp = The_list_of_all_edges->data[j];
                The_list_of_all_edges->data[j] = The_list_of_all_edges->data[j + 1];
                The_list_of_all_edges->data[j + 1] = temp;
            }
        }
    }
}
static int find(int parent_track[], int vertexno)
{
    return (parent_track[vertexno]); // I have used this function finds the absoulte root
}
static void applyUnion(int n, int parent_track[], int c1, int c2)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (parent_track[i] == c2) // I have used this function to apply the union operaion
        {
            parent_track[i] = c1;
            // break;
        }
    }
}

static void Kruskals_Algorithim(int n, pair_t edges[n - 1], const connection_t connections[n][n], edge_list *The_list_of_all_edges, edge_list *The_span_list)
{
    int parent_track[n];
    int i;
    int absolute_root1;
    int absolute_root2;
    bubble_sort(The_list_of_all_edges); // This function sorts the lis of edges based on time using bubble sort

    for (i = 0; i < n; i++)
    {
        parent_track[i] = i; // This array keeps track of which set a edge in before we add a new one
    }

    The_span_list->n = 0; // this is the span list

    for (i = 0; i < The_list_of_all_edges->n; i++) // selects each edge from the sorted array
    {
        absolute_root1 = find(parent_track, The_list_of_all_edges->data[i].u); // finds the absoulute node of u
        absolute_root2 = find(parent_track, The_list_of_all_edges->data[i].v); //

        if (absolute_root1 != absolute_root2) // check if the node is not in the set(checks if a loop is created beofre addig)
        {
            The_span_list->data[The_span_list->n] = The_list_of_all_edges->data[i]; // adding the edge to The_span_list
            The_span_list->n = The_span_list->n + 1;                                // span list increases
            applyUnion(n, parent_track, absolute_root1, absolute_root2);            // actually attaching the two sets
        }
    }
}
static int min_two(int x, int y)
{
    if (x < y)
    {
        return x;
    }
    else
    {
        return y;
    }
}

static int get_all_trips(int n,
                         const connection_t connections[n][n],
                         int visited[],
                         int start,
                         int currPos,
                         int exempted,
                         int visit_count,
                         int cost,
                         trip all_trips[n],
                         int *size)
{

    visited[currPos] = visit_count;
    visit_count++;

    if (visit_count == (n - 1))
    {
        // (n-2) is the maximum we are allowed to visit, excluding
        // start airport and the exempted airport.

        if (connections[currPos][start].distance == INT_MAX)
        {
            // After having visited (n-2) airports, the airport
            // where we are at now (currPos) does not have a connection to
            // the start airport.

            return 0; // Failure
        }
        else
        {
            /*stores the total cost*/
            all_trips[*size].cost = (cost + connections[currPos][start].distance);

            extract_trip_seq(n - 1, visited, all_trips[*size].sequence_of_visits);

            *size = *size + 1;
            return 1; // Success
        }
    } // if we have visited (n-2) airports so far.

    // BACKTRACKING STEP
    // Loop to traverse the adjacency list
    // of currPos node and increasing the visit_count
    // by 1 and cost by graph[currPos][i] value
    for (int i = 0; i < n; i++)
    {
        /* code */
        if ((i != exempted) && (i != start))
        {
            if (visited[i] == 0 && connections[currPos][i].distance != INT_MAX) // work
            {
                // Mark as visited
                visited[i] = visit_count;
                int ret_val;

                ret_val = get_all_trips(n,
                                        connections,
                                        visited,
                                        start,
                                        i,
                                        exempted,
                                        visit_count,
                                        cost + connections[currPos][i].distance,
                                        all_trips, size);
                // Mark ith node as unvisited
                visited[i] = 0;
            }
        }
    } // loop through adjacency list.

    return 0; // Failure
}

static int min(int n, trip all_trips[n])
{
    int min = 0;
    for (int i = 0; i < n; i++)
    {
        if (all_trips[i].cost < all_trips[min].cost)
        {
            min = i;
        }
    }
    return min;
}
static int find_min_time(int n,int v,int dst,int curCost,int minCost,const connection_t connections[n][n],int visited[n])
{
    // int res = INT_MAX;
    visited[v] = 1;
    if (v == dst)
    {
        if (curCost < minCost)
        {
            minCost = curCost;
        }
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            if (!visited[i] && connections[v][i].time != INT_MAX && connections[v][i].time != 0)
            {
                minCost = find_min_time(n, i, dst, curCost + connections[v][i].time, minCost, connections, visited);
            }
        }
    }
    visited[v] = 0;
    return minCost;
}

// YOUR SOLUTIONS BELOW

int q1(int n, const connection_t connections[n][n])
{ ////so for this question, I have decided to use a modified version of kosaraju's algorithim which also utilizes dfs traversal
    int airports[n];
    for (int i = 0; i < n; i++)
    {
        airports[i] = 0; // extracting the list of aiports this will also behave like a visited arrray for DFS
    }
    int stack[n];
    int stack1[n];
    int top = 0;
    stack[n] = '\0';

    connection_t rev_graph[n][n];
    reverse_graph(n, connections, rev_graph);      // works
    dfs(0, n, connections, airports, stack, &top); // wow it works
    for (int i = 0; i < n; i++)
    {
        airports[i] = 0; // reseting the airports1 array
    }
    top = 0;                                      // reseting the top
    dfs(0, n, rev_graph, airports, stack1, &top); // now doing a dfs travel from the same start point and storing the path in stack1
    int flag;
    if (count_ones(n, stack) == count_ones(n, stack1)) // checking the total nnumber of visited nodes in both the graphs and comparing them
    {
        return 1;
    }
    else
    {
        flag = 0;
    }

    // It worked yaaaaaaaaaaaaaay

    return flag;
}

int q2(const airport_t *src, const airport_t *dest, int n, int k, const connection_t connections[n][n])
{
    // so basically we need to find all paths
    int visited[n];
    for (int i = 0; i < n; i++)
    {
        visited[i] = 0; // initializing the visited array
        // which be useful the upcomming recursive path_find function
    }
    int flag;
    flag = path_find(n, dest->num_id, connections, k, src->num_id, visited, 0); // recursive call begins
    if (flag == 1)
    {
        return 1;
    }

    return 0;
}

int q3(const airport_t *src, int n, const connection_t connections[n][n])
{
    // Finding a loop in a graph
    // for that we can re-use path_find() function!
    int adjacent[n];
    int end = 0;
    for (int j = 0; j < n; j++)
    {
        if ((connections[src->num_id][j].distance != INT_MAX) && (connections[src->num_id][j].time != INT_MAX))
        {
            adjacent[end] = j; /* we are collecting all the airports adjacent  to  the source's
                                */
            end = end + 1;
        }
    }
    int visited[n];

    for (int i = 0; i < end; i++)
    {
        for (int i = 0; i < n; i++)
        {
            visited[i] = 0; // initializing the visited array
        }
        if (path_find(n, src->num_id, connections, n, adjacent[0], visited, 0))
        {
            return 1; // path found
        }
    }
    // if we can't find a path from anyone of  the adjacent-to-src nodes
    // to src we come here
    return 0;
}

void q4(int n, int (*predicate_func)(const airport_t *, const airport_t *),
        airport_t airport_list[n])
{

    quickSort(n, airport_list, predicate_func, 0, n - 1); // IT WORKED!!!!!!!!
    // yaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaayy!
}

pair_t q5(int n, airport_t airports[n])
{
    char common_prefix[n];
    int len = 0;
    int j;
    pair_t ans = {-1, -1};

    for (int i = 0; i < n; i++)
    {
        if (len < definently_not_strlen(airports[i].airport_name)) // skiping a few unnecessary comparissons
        {
            j = i + 1;
            while (j < n)
            {

                int k;
                int flag;
                for (k = 0, flag = 1; flag &&
                                        ((airports[i].airport_name[k] != '\0') &&
                                         (airports[j].airport_name[k] != '\0')); )
                { // the loop && following if block will compare two strings
                    if (airports[i].airport_name[k] != airports[j].airport_name[k])
                    {
                        flag = 0;
                    }
                    else
                    {
                        k += 1;
                    }
                }
                if (k > (len + 1))
                {
                    // common_prefix[k] = airports[i].airport_name[k];//src=airports[i].airport_name[k]; dest=commonprefix,n=k-1
                    definently_not_strncpy(common_prefix, airports[i].airport_name, k);

                    len = k;
                    ans.first = i;
                    ans.second = j;
                }

                j = j + 1;
            }
            // j=i+1;
        }
    }

    return ans;
}

int q6(int n, int amount, const int entry_fee[n])
{

    // Driver code

    int ans = binarySearchCount(n, entry_fee, amount);
    // huh it actually correctly on first go nice!
    return ans;
}

void q7(int n, const char *pat, int contains[n], const airport_t airports[n])
{
    int where; // this variable is the index of where the pattern is found
    for (int i = 0; i < n; i++)
    {

        where = horspool(airports[i].airport_name, pat); // while this one does the actual checking
        if (where >= 0)
        {
            contains[i] = 1;
        }
        else
        {
            contains[i] = 0;
        }
    }
}

int q8(int n, int trip_order[n - 1], const connection_t connections[n][n])
{
    int visited[n];

    trip all_trips[5040]; // factorial(7)
    int size = 0;

    if (n <= 2)
    {
        // We need a minimum of 3 airports, so that one can be the
        // start airport, one can be the exempted airport, and the
        // remaining one can be the visited airport.

        return -1;
    }

    // trip_order needs to be passed
    for (int i = 0; i < n; i++) // the start-point  loop
    {
        for (int k = 0; k < n; k++) // the exemption loop
        {
            for (int j = 0; j < n; j++) // initializing visited array
            {
                visited[j] = 0;
            }
            if ((k != i) && (connections[i][k].distance != INT_MAX))
            {
                visited[k] = -1; // mark k as exempted
                get_all_trips(n, connections, visited, i, i, k, 0, 0, all_trips, &size);
            }
        } // the k loop
    }     // the i loop

    if (size > 0) // If at least one full-coverage cyclic trips found
    {
        int min_index = min(n, all_trips);
        for (int i = 0; i < n - 1; i++)
        {
            trip_order[i] = all_trips[min_index].sequence_of_visits[i];
        }
        int min_cost = all_trips[min_index].cost;
        return min_cost;
    }
    else
    { // not changing trip_order
        return -1;
    }
}

int q9(int n, pair_t edges[n - 1], const connection_t connections[n][n])
{
    edge_list The_list_of_all_edges;
    edge_list The_span_list;
    The_list_of_all_edges.n = 0;
    The_span_list.n = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (connections[i][j].time != 0 && connections[i][j].time != INT_MAX)
            {
                The_list_of_all_edges.data[The_list_of_all_edges.n].u = i;
                The_list_of_all_edges.data[The_list_of_all_edges.n].v = j;
                The_list_of_all_edges.data[The_list_of_all_edges.n].w = connections[i][j]; // collects all the info of all the edges
                The_list_of_all_edges.n++;
            }
        }
    }
    Kruskals_Algorithim(n, edges, connections, &The_list_of_all_edges, &The_span_list);
    int count = 0;
    for (int i = 0; i < n - 1; i++)
    {
        count = count + The_span_list.data[i].w.time;
        edges[i].first = The_span_list.data[i].u;
        edges[i].second = The_span_list.data[i].v;
    }

    return count;
}

void q10(int n, int k, const airport_t *src,
         const connection_t connections[n][n], const int destinations[k],
         int costs[k])
{
    int visited[n];
    int source = src->num_id;

    for (int i = 0;i < k;++i)
    {
        for(int i=0;i<n;i++)
        {
            visited[i]=0;
        }
        costs[i] = find_min_time(n,source,destinations[i],0,INT_MAX,connections,visited);
    }
}

// END
