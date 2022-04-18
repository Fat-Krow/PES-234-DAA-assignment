#include "header.h"
#include<limits.h>
static void reverse_graph (int n,const connection_t connection[n][n],connection_t rev[n][n]);//returns a transposed digraph i.e reverses the edges in the digraph
static void dfs (int i,int n,const connection_t connection[n][n],int airports[n],int stack[n],int* top);
int count_ones(int n,int arr[n]);
int path_find(int n,int  dest ,const connection_t connection[n][n] ,int k,int sour,int *visited,int curr_path_length);


// ANY STATIC FUNCTIONS ARE UP HERE
static void  reverse_graph(int n, const connection_t connection[n][n],connection_t rev[n][n])
{   //connection_t copy;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            rev[i][j]=connection[j][i];

        }

    }


}
static void dfs (int i,int n,const connection_t connection[n][n],int airports[n],int stack[n],int* top)
{
    airports[i]=1;

    for(int j=0; j<n ;j++)
    {   //edge=connection[i][j];
        if ( (connection[i][j].distance!=INT_MAX)&&(connection[i][j].time!=INT_MAX))//error
        {   if (airports[j]==0)
                {dfs(j,n,connection,airports,stack,top);}
        }
    }

    stack[*top]=i;
    *top=*top+1;


}
int count_ones(int n,int arr[n])
{   int count=0;
    for (int i=0;i<n;i++)
    {
        if(arr[i]==1)
        {
            count=count+1;
        }
    }
    return count;
}
int path_find(int n,int  dest ,const connection_t connection[n][n] ,int k,int src,int *visited,int curr_path_length)
{
    if(src==dest)
    {
        return 1;
    }

    if(k<curr_path_length)
    {
        return 0;//test
    }
    visited[src]=1;

        for(int j=0;j<n;j++)
        {
            if(visited[j]!=1 && connection[src][j].distance!=INT_MAX && connection[src][j].time!=INT_MAX )
            {
               int flag= path_find(n,dest,connection,k,j,visited,curr_path_length+1);

               if (flag==1)
                {
                        return 1;
                }


            }

        }
        return 0;//from this node there is no unvisted node

}

// YOUR SOLUTIONS BELOW

int q1(int n, const connection_t connections[n][n])
{   ////so for this question, I have decided to use a modified version of kosaraju's algorithim which also utilizes dfs traversal
    int airports[n];
    for (int i = 0; i < n; i++)
    {
        airports[i]=0;//extracting the list of aiports this will also behave like a visited arrray for DFS
    }
    int stack[n];
    int stack1[n];
    int top=0;
    stack[n]='\0';

    connection_t rev_graph [n][n];
    reverse_graph(n,connections,rev_graph);//works
    dfs(0,n,connections,airports,stack,&top);//wow it works
    for (int i = 0; i < n; i++)
    {
        airports[i]=0;//reseting the airports1 array
    }
    top=0;//reseting the top
    dfs(0,n,rev_graph,airports,stack1,&top);//now doing a dfs travel from the same start point and storing the path in stack1
    int flag;
    if(count_ones(n,stack)==count_ones(n,stack1))//checking the total nnumber of visited nodes in both the graphs and comparing them
    {
        return 1;
    }
    else
    {
        flag=0;
    }

//It worked yaaaaaaaaaaaaaay

    return flag;
}

int q2(const airport_t *src, const airport_t *dest, int n, int k,const connection_t connections[n][n])
{
    //so basically we need to find all paths
    int visited[n];
    for (int i=0;i<n;i++)
    {
        visited[i]=0;
    }
    int flag;
    flag=path_find(n,dest->num_id,connections,k,src->num_id,visited,0);
    if (flag==1)
    {
        return 1;

    }


    return 0;
}

int q3(const airport_t *src, int n, const connection_t connections[n][n])
{



    return 0;
}

void q4(int n, int (*predicate_func)(const airport_t *, const airport_t *),
        airport_t airport_list[n])
{
}

pair_t q5(int n, airport_t airports[n])
{
    pair_t ans = {-1, -1};
    return ans;
}

int q6(int n, int amount, const int entry_fee[n])
{
    return 0;
}

void q7(int n, const char *pat, int contains[n], const airport_t airports[n])
{
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
