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
static int horspool(const char src[],const char p[]);
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
        return 1;
    }

    if (k < curr_path_length)
    {
        return 0; // test
    }
    visited[src] = 1;

    for (int j = 0; j < n; j++)
    {
        if (visited[j] != 1 && connection[src][j].distance != INT_MAX && connection[src][j].time != INT_MAX)
        {
            int flag = path_find(n, dest, connection, k, j, visited, curr_path_length + 1);

            if (flag == 1)
            {
                return 1;
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

    // select the rightmost element as pivot
    airport_t pivot = array[high];

    // pointer for greater element
    int i = (low - 1);

    // traverse each element of the array
    // compare them with the pivot
    for (int j = low; j < high; j++)
    {
        if (predicate_func(&array[j], &pivot))
        {

            // if element smaller than pivot is found
            // swap it with the greater element pointed by i
            i++;

            // swap element at i with element at j
            swap(&array[i], &array[j]);
        }
    }

    // swap the pivot element with the greater element at i
    swap(&array[i + 1], &array[high]);

    // return the partition point
    return (i + 1);
}
static void quickSort(int n, airport_t array[n], int (*predicate_func)(const airport_t *, const airport_t *), int low, int high)
{

    if ((low < high) && predicate_func(&array[low], &array[high])) // this is for ascending
    {

        // find the pivot element such that n-1
        // elements smaller than pivot are on left of pivot
        // elements greater than pivot are on right of pivot
        int pi = partition(n, array, predicate_func, low, high);

        // recursive call on the left of pivot

        quickSort(n, array, predicate_func, low, pi - 1);

        // recursive call on the right of pivot
        quickSort(n, array, predicate_func, pi + 1, high);
    }
    else if ((low < high) && predicate_func(&array[high], &array[low])) // this for descending patterns
    {

        // find the pivot element such that n-1
        // elements smaller than pivot are on left of pivot
        // elements greater than pivot are on right of pivot
        int pi = partition(n, array, predicate_func, low, high);

        // recursive call on the left of pivot

        quickSort(n, array, predicate_func, low, pi - 1);

        // recursive call on the right of pivot
        quickSort(n, array, predicate_func, pi + 1, high);
    }
}
static void definently_not_strncpy(char *dest, const char *src, int k)
{
    int i;
    for (i = 0; ((i < k) && (src[i] != '\0')); i++)
    {
        dest[i] = src[i];
    }
    dest[i] = '\0';
}
static int definently_not_strlen(const char *string)
{
    unsigned int count = 0;
    while(*string!='\0')
    {
        count++;
        string++;
    }
    return count;
}
// A binary search function to return
// the number of elements less than
// or equal to the given key
static int binarySearchCount(int n, const int arr[n], int key)
{

    int left = 0;
    int right = n - 1;

    int count = 0;
    int mid;
    while (left <= right)
    {
        mid = (right + left) / 2;

        // Check if middle element is
        // less than or equal to key
        if (arr[mid] <= key)
        {
            // At least (mid + 1) elements are there
            // whose values are less than
            // or equal to key
            count = mid + 1;
            left = mid + 1;
        }

        // If key is smaller, ignore right half
        else
        {
            right = mid - 1;
        }
    }

    return count;
}
static void shifttable(int shift_table[256],const char pattern[])
{ /*this function will create the shift table */
    int i, j, m;
    m = definently_not_strlen(pattern); // m=length of pattern
    for (i = 0; i < 256; i++)
    {
        shift_table[i] = m;
    }
    for (j = 0; j < m - 1; j++)
    {
        shift_table[(int)pattern[j]] = m - 1 - j;
    }
}
static int horspool(const char text[],const char pattern[])
{   int  shift_table[256];
    shifttable(shift_table,pattern);
    int i,k, m, n;
    n = definently_not_strlen(text);//length of the string IN WHICH we need to search
    m = definently_not_strlen(pattern);//length of the pattern we need to search for
    i = m - 1;//because we start from last char
    //int where;
    //shifttable(shift_table,pattern);//this function prepares the shift table
    while (i < n)
    {
        k = 0;//first char of pattern string
        
        while ((k < m) && (pattern[m - 1 - k] == text[i-k]))//i=2-1
        {
            k=k+1;//the real horspool logic starts here
        }
        if (k == m)
        {
             
            return (i - m + 1);//returns where it is found
        }
        else
        {
                i = i + shift_table[(int)text[i]];
        }
    }
    return -1;
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
        visited[i] = 0;
    }
    int flag;
    flag = path_find(n, dest->num_id, connections, k, src->num_id, visited, 0);
    if (flag == 1)
    {
        return 1;
    }

    return 0;
}

int q3(const airport_t *src, int n, const connection_t connections[n][n])
{
    // Finding a loop in a graph
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
            visited[i] = 0;
        }
        if (path_find(n, src->num_id, connections, n, adjacent[0], visited, 0))
        {
            return 1;
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
                for (k = 0; (airports[i].airport_name[k] != '\0' && airports[j].airport_name[k] != '\0'); k++)
                { // the loop && following if block will compare two strings
                    if (airports[i].airport_name[k] != airports[j].airport_name[k])
                    {
                        break;
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

    return ans;
}

void q7(int n, const char *pat, int contains[n], const airport_t airports[n])
{
    int where;
    for (int i = 0; i < n; i++)
    {

        where=horspool(airports[i].airport_name,pat);//while this one does the actual checking
        if(where>=0)
        {
            contains[i]=1;
        }
        else
        {
            contains[i]=0;
        }

    }

}

int q8(int n, int trip_order[n - 1], const connection_t connections[n][n])
{
    return 0;
}

int q9(int n, pair_t edges[n - 1], const connection_t connections[n][n])
{
    return 0;
}

void q10(int n, int k, const airport_t *src,
         const connection_t connections[n][n], const int destinations[k],
         int costs[k])
{
}

// END
