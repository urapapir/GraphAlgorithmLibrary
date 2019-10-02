
#pragma once
#include <vector>
#include <iostream>
#include<fstream>
#include <iomanip>
#include <conio.h>
#include <algorithm>
#include <map>
#include <functional>
#include <math.h>
#include<stack>
#include <iterator> 
using namespace std;

struct e {
	int c, a;
};
typedef vector<vector<int>> Graph;
typedef vector < vector < e > > Network;
typedef vector<vector<int>> Graph;
typedef vector<pair<int, pair<int, int>>> Graph_edges; //(вес, (вершина 1, вершина 2))
const int inf = 1000000000;

namespace GraphAlgorithmLibrary
{
	int min(int a, int b) {
		if (a == 0) return b;
		if (b == 0) return a;
		return (a < b) ? a : b;
	}

class GraphAlgorithmLibrary
{

		int iter = 0;
		struct e {
			int c, a;
		};
		typedef vector<vector<int>> Graph;
		typedef vector < vector < e > > Network;
		typedef vector<vector<int>> Graph;
		typedef vector<pair<int, pair<int, int>>> Graph_edges; //(вес, (вершина 1, вершина 2))

		Graph Graph;
		Network Network;
		ofstream out();

		GraphAlgorithmLibrary() {
			Graph graph(0, vector<int>(0));
			ofstream out("output");
		}

		Graph GenerateGraph(int n, double p = 0.5, int min_weight = -100, int max_weight = 100, double p_weight = 0.5);
		Network GenerateNetwork(const Graph& graph);
		Graph ReadGraph(const char* name);

		int determinant(Graph graph); 
		vector <int> ddeg(Graph general);
		Graph To_undirected(const Graph& graph);
		Graph_edges Matrix_to_Edges(const Graph& graph);
		bool noway(const Graph& matr, int p, int tv, int in);//проверяет есть ли пути до заданных вершин
		void Reachable(const Graph& graph, int a, int b, ostream& out = cout, bool Print = false);
		Graph Hamiltonian(const Graph& graph);
		vector<int> GenerateBinomial(int m, int n, double p = 0.5);


		void Floyd(const Graph& graph, int a, int b);
		bool Shimbell_new(const Graph& graph, int count, bool Print = true, int a = 0, int b = 0, ostream & out = cout);
		vector<int> Dijkstra_new(const Graph& graph, int a, int b, bool print = true);
		vector<int> Bellman_Ford(const Graph& graph, int a, int b, bool print = true);
		void Prim_algorithm(const Graph& graph);
		void Kruskal_algorithm(const Graph& graph);
		Graph Prim(const Graph& graph);		
		void Kirchhoff(Graph graph);
		vector<int> PruferCode(Graph graph);
		vector<int> Eiler(vector<vector<int>> matr);
		int Salesman(const Graph& graph, vector<int>& cycle, int v = 0);//поиск минимальной длины гамильтонова цикла
		void Hamilton(const Graph& graph);
		int Fulkerson(const Network& network, int s, int t);
        int FindMinCostFlow(const Network& network, int s, int t, int flow);


		void InfToZeros(Graph& graph);
		void ZerosToInf(Graph& graph);
		bool isEilerGraph(const vector<vector<int>>& graph);
		bool isNull(const Graph& graph);
		bool Svyazn(const Graph& matr, int tv, int in);
		
		
		friend ostream& operator <<(ostream& out, const Graph& graph);
		friend ostream& operator <<(ostream& out, const Network& graph);
		


	};

}

