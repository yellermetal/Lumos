
/*
*	HopcroftKarp implementation
*	Language: C
*	Author: Ariel Livshits 2018
*/

/* Includes */
#include "HopcroftKarp.h"

/* Typedefs and Defines */
#define NIL -1

/* Declarations */
void fill(int_ptr arr, int size, int val);
void fill(bool_ptr arr, int size, bool val);
void freeHopcroftKarp(HopcroftKarp_ptr graph);
void addEdge(int u, int v, HopcroftKarp_ptr graph);
bool dfs(int u1, HopcroftKarp_ptr graph);
void bfs(HopcroftKarp_ptr graph);
int maxMatching(HopcroftKarp_ptr graph);


/* ------------------------------------------------- */

void HopcroftKarp(lms_mat_ptr mat_ptr, int_ptr* out_matching, int_ptr out_maxcardinality)
{
	if (!mat_ptr) return;

	HopcroftKarp_ptr graph = (HopcroftKarp_ptr)malloc(sizeof(HopcroftKarp_t));
	graph->n = mat_ptr->dim;
	int size = (graph->n)*(graph->n);
	size_t int_size = sizeof(int)*size;
	size_t bool_size = sizeof(bool)*size;

	graph->last = (int_ptr)malloc(int_size); graph->prev = (int_ptr)malloc(int_size); graph->head = (int_ptr)malloc(int_size);
	graph->matching = (int_ptr)malloc(int_size); graph->dist = (int_ptr)malloc(int_size); graph->Q = (int_ptr)malloc(int_size);
	graph->used = (bool_ptr)malloc(bool_size); graph->vis = (bool_ptr)malloc(bool_size);

	graph->edges = 0;
	fill(graph->last, graph->n, NIL);

	int_ptr nz_elems = getNonzeroElemsIndices(mat_ptr);

	int left_vertex;
	int right_vertex;
	for (int idx = 0; idx < mat_ptr->nz; idx++) {
		left_vertex = mat_ptr->elems[nz_elems[idx]].row;
		right_vertex = mat_ptr->elems[nz_elems[idx]].col;
		addEdge(left_vertex, right_vertex, graph);

	}
	free(nz_elems);
	
	*out_maxcardinality = maxMatching(graph);

	if (*out_maxcardinality == 0 || out_matching == NULL) {
		if (out_matching != NULL) out_matching = NULL;
		free(graph->matching);
		freeHopcroftKarp(graph);
		return;
	}

	*out_matching = graph->matching;

	freeHopcroftKarp(graph);

}

void freeHopcroftKarp(HopcroftKarp_ptr graph)
{
	free(graph->last); 
	free(graph->prev); 
	free(graph->head);
	free(graph->dist);
	free(graph->Q); 
	free(graph->used); 
	free(graph->vis);
	free(graph);
}
void addEdge(int u, int v, HopcroftKarp_ptr graph)
{
	graph->head[graph->edges] = v;
	graph->prev[graph->edges] = graph->last[u];
	graph->last[u] = graph->edges++;
}
bool dfs(int u1, HopcroftKarp_ptr graph)
{
	graph->vis[u1] = true;
	for (int e = graph->last[u1]; e >= 0; e = graph->prev[e]) {
		int v = graph->head[e];
		int u2 = graph->matching[v];
		if (u2 < 0 || !graph->vis[u2] && graph->dist[u2] == graph->dist[u1] + 1 && dfs(u2, graph)) {
			graph->matching[v] = u1;
			graph->used[u1] = true;
			return true;
		}
	}
	return false;
}
void bfs(HopcroftKarp_ptr graph)
{
	fill(graph->dist, graph->n, NIL);
	int sizeQ = 0;
	for (int u = 0; u < graph->n; ++u) {
		if (!graph->used[u]) {
			graph->Q[sizeQ++] = u;
			graph->dist[u] = 0;
		}
	}
	for (int i = 0; i < sizeQ; i++) {
		int u1 = graph->Q[i];
		for (int e = graph->last[u1]; e >= 0; e = graph->prev[e]) {
			int u2 = graph->matching[graph->head[e]];
			if (u2 >= 0 && graph->dist[u2] < 0) {
				graph->dist[u2] = graph->dist[u1] + 1;
				graph->Q[sizeQ++] = u2;
			}
		}
	}
}
int maxMatching(HopcroftKarp_ptr graph)
{
	fill(graph->used, graph->n, false);
	fill(graph->matching, graph->n, NIL);
	for (int res = 0;;) {
		bfs(graph);
		fill(graph->vis, graph->n, false);
		int f = 0;
		for (int u = 0; u < graph->n; ++u)
			if (!graph->used[u] && dfs(u, graph))
				++f;
		if (!f)
			return res;
		res += f;
	}
}
void fill(int_ptr arr, int size, int val) {
	for (int iter = 0; iter < size; iter++)
		arr[iter] = val;
}
void fill(bool_ptr arr, int size, bool val) {
	for (int iter = 0; iter < size; iter++)
		arr[iter] = val;
}
