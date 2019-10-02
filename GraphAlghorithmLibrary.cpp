#include "stdafx.h"
#include "GraphAlghorithmLibrary.h"
//const int inf = 1000000000;
namespace GraphAlgorithmLibrary
{

vector<int> Chetnost(Graph matr) {
	int size = matr.size();
	vector<int> Chetn(size, 0);
	for (int i = 0; i < matr.size(); i++)
		for (int j = 0; j < matr.size(); j++)
			if (matr[i][j])
				Chetn[i]++;
	return Chetn;
}

bool noway(const Graph& matr, int p, int tv, int in) {//��������� ���� �� ���� �� �������� ������
	int n = matr.size();
	for (int i = 0; i < n; i++)
		if (matr[p][i] || p == in || p == tv) return true;
	return false;
}//��������� ���� �� ���� �� �������� ������
int determinant(Graph graph) {
	int i, j;
	int N = graph.size();
	int determ = 0;
	Graph graph1 = graph;
	for (int i = 0; i < graph.size(); i++) {
		for (int j = 0; j < graph.size(); j++) {
			if (graph1[i][j] == inf) { graph1[i][j] = 0; }
		}
	}
	Graph matr1(0);

	if (N == 1)
	{
		determ = graph1[0][0];
	}
	else
	{
		matr1.resize(N - 1);

		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N - 1; j++)
			{
				if (j < i)
				{
					matr1[j] = graph1[j];
				}
				else
				{
					matr1[j] = graph1[j + 1];
				}
			}
			determ += pow(-1.0, (i + j)) * determinant(matr1) * graph1[i][N - 1];
		}
	}

	return determ;
}

vector <int> ddeg(Graph general) {//���������� ������� ������
	vector <int> deg;
	for (int i = 0; i < general.size(); i++) {
		deg.push_back(0);
		for (int j = 0; j < general.size(); j++) {
			//deg[i] += (general[i][j]) ? 1 : 0;
			if (general[i][j] != inf && general[i][j] != 0) {
				deg[i] += 1;
			}
		}
	}
	return deg;
}

ostream& operator <<(ostream& out, const Graph& graph) {
	out << "     ";
	for (int i = 0; i < graph.size(); i++) {
		out << setw(3) << i + 1 << " ";
	}
	out << endl;
	for (int i = 0; i < graph.size(); i++) {
		out << setw(3) << i + 1 << "  ";
		for (int j = 0; j < graph.size(); j++) {
			if (graph[i][j] < inf) {
				out << setw(3) << (graph[i][j]) << " ";
			}
			else {
				out << setw(3) << "*" << " ";
			}
		}
		out << endl;
	}
	return out;
}

ostream& operator <<(ostream& out, const Network& graph) {
	cout << "������� ����������  ������������: " << endl;
	out << "     ";
	for (int i = 0; i < graph.size(); i++) {
		out << setw(3) << i + 1 << " ";
	}
	out << endl;
	for (int i = 0; i < graph.size(); i++) {
		out << setw(3) << i + 1 << "  ";
		for (int j = 0; j < graph.size(); j++) {
			out << setw(3) << graph[i][j].c << " ";
		}
		out << endl;
	}
	cout << endl << "������� ����������: " << endl;
	out << "     ";
	for (int i = 0; i < graph.size(); i++) {
		out << setw(3) << i + 1 << " ";
	}
	out << endl;
	for (int i = 0; i < graph.size(); i++) {
		out << setw(3) << i + 1 << "  ";
		for (int j = 0; j < graph.size(); j++) {
			out << setw(3) << graph[i][j].a << " ";
		}
		out << endl;
	}
	return out;
}

Graph ReadGraph(const char* name) {
	ifstream in(name);
	int n;
	in >> n;
	int tmp;
	for (int i = 1; i <= n; i++) {

		in >> tmp;
	}
	Graph graph(n, vector<int>(n));
	for (int i = 0; i < n; i++) {
		in >> tmp;
		for (int j = 0; j < n; j++) {
			in >> graph[i][j];
		}
	}
	return graph;
}

Graph_edges Matrix_to_Edges(const Graph& graph) {//���������� ������� �����
	Graph_edges edges;
	int n = graph.size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (graph[i][j] < inf) {
				edges.push_back(make_pair(graph[i][j], make_pair(i, j)));
			}
		}
	}
	return edges;
}

vector<int> Bellman_Ford(const Graph& graph, int a, int b, bool print = true) {
	int n = graph.size();
	if (a < 1 || a > n || b < 1 || b > n) {
		cout << "Graph vertex out of range." << endl;
		return vector<int>();
	}
	a--; b--;
	Graph_edges edges = Matrix_to_Edges(graph);
	int m = edges.size();
	vector<int> Distances(n, inf), Prev(n, -1);
	Distances[a] = 0;
	Prev[b] = a;
	bool changed;
	do {
		changed = false;
		for (int i = 0; i < m; i++) {
			//edges[i] - ����� �� u � v ���� weight
			int u = edges[i].second.first, v = edges[i].second.second, weight = edges[i].first;//u - 1-�� �������, v - ������
			if (Distances[u] < inf && Distances[u] + weight < Distances[v]) {
				Distances[v] = Distances[u] + weight;//���������� ����� (u,v)
				Prev[v] = u;
				changed = true;
			}
		}
	} while (changed);
	vector<int> Path;
	for (int i = b; i != -1; i = Prev[i]) {
		Path.push_back(i);
	}
	reverse(Path.begin(), Path.end());
	if (print) {
		if (Distances[b] == inf) {
			cout << "��� ���� �� " << a + 1 << " � " << b + 1 << ":\n";
		}
		else {
			cout << "����������� ���� �� " << a + 1 << " �� " << b + 1 << ":\n";
			for (int i = 0; i < Path.size(); i++) {
				cout << Path[i] + 1 << " ";
			}
			cout << endl << "�����: " << Distances[b] << ".\n";
		}
	}
	return Path;
}

void Floyd(const Graph& graph, int a, int b) {
	int n = graph.size();
	if (a < 1 || a > n || b < 1 || b > n) {
		cout << "Graph vertex out of range." << endl;
		return;
	}
	a--; b--;
	Graph Distances = graph, Next(n, vector<int>(n, -1));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			//if (!Distances[i][j]){
			//	Distances[i][j] = inf;
			//}
			Next[i][j] = j;
		}
	}
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (Distances[i][k] < inf && Distances[k][j] < inf && Distances[i][k] + Distances[k][j] < Distances[i][j]) {
					Distances[i][j] = Distances[i][k] + Distances[k][j];
					Next[i][j] = Next[i][k];
				}
			}
		}
	}
	if (Distances[a][b] == inf) {
		cout << "��� ���� �� " << a + 1 << " � " << b + 1 << ":\n";
	}
	else {
		cout << "����������� ���� �� " << a + 1 << " �� " << b + 1 << ":\n";
		for (int i = a;; i = Next[i][b]) {
			cout << i + 1 << " ";
			if (i == Next[i][b]) {
				break;
			}
		}
		cout << endl << "�����: " << Distances[a][b] << ".\n";
	}
}

vector<int> GenerateBinomial(int m, int n, double p = 0.5) {//���������� �� ������������ �������
	vector<double> C(n + 1), C0(n + 1);//������������ �����������
	vector<double> P(n + 1), Q(n + 1);
	C0[0] = 1;//n=0
	C[0] = 1;//n=1
	P[0] = 1;
	Q[0] = 1;
	for (int i = 1; i <= n; i++) {
		C0.swap(C);
		for (int k = 1; k <= i; k++) {
			C[k] = C0[k - 1] + C0[k];
		}
		P[i] = P[i - 1] * p;
		Q[i] = Q[i - 1] * (1 - p);
	}
	vector<double> Probabilities(n);
	for (int k = 0; k < n; k++) {
		Probabilities[k] = C[k] * P[k] * Q[n - k];//������� �����������
	}
	vector<int> result(m);
	for (int i = 0; i < m; i++) {
		double a = double(rand()) / RAND_MAX;
		int k;
		for (k = 0; k < n && a >= Probabilities[k]; k++) {
			a -= Probabilities[k];
		}
		result[i] = k;
	}
	return result;
}

Graph GenerateGraph(int n, double p = 0.5, int min_weight = -100, int max_weight = 100, double p_weight = 0.5)
{
	if (n < 0) {
		cout << "Nepravil'no vvedeno" << endl;
		n = 0;
	}
	Graph graph(n, vector<int>(n, inf));//graph
	if (n == 0) {
		return graph;
	}
	vector<int> Weights = GenerateBinomial(n - 1, max_weight - min_weight, p_weight);
	for (int i = 0; i < n - 1; i++) {
		graph[i][i + 1] = Weights[i] + min_weight;
	}
	if (n <= 2) {
		return graph;
	}
	vector<int> Tmp = GenerateBinomial(n, n - 3, p);
	Tmp[0]++;
	Tmp[n - 1]++;
	int M = 0;
	multimap<int, int, greater<int>> Degrees;//�������/����� �������
	for (int i = 0; i < n; i++) {
		if (Tmp[i] != 0)
			Degrees.insert(make_pair(Tmp[i], i));
	}
	for (int i = 0; i < n; i++) {
		M += Tmp[i];
	}
	M /= 2;
	Weights = GenerateBinomial(M, max_weight - min_weight, p_weight);
	int index = 0;
	while (!Degrees.empty()) {
		multimap<int, int, greater<int>>::iterator i = Degrees.begin(), i0 = i;
		int v = i->second;
		int d = i->first;//���������� ���������� ��������� ������� � �����
		for (++i; i != Degrees.end() && d > 0; ++i) {//����. ����� ������� � ������ ���-���
			int u = i->second;
			if (graph[v][u] == inf && graph[u][v] == inf) {
				if (v < u) {
					graph[v][u] = Weights[index] + min_weight;
				}
				else {
					graph[u][v] = Weights[index] + min_weight;
				}
				index++;
				d--;
				int d1 = i->first - 1;
				Degrees.erase(i);
				if (d1 > 0) {
					Degrees.insert(make_pair(d1, u));
				}
				i = i0;
			}
			else {
				i0 = i;
			}
		}
		Degrees.erase(Degrees.begin());
	}
	for (int i = 0; i < graph.size(); i++) {
		for (int j = 0; j < graph.size(); j++) {
			if (graph[i][j] == 0) {
				graph[i][j] = 3;
			}
		}
	}
	return graph;
}

bool Shimbell(const Graph& graph, int count, bool Print = true, int a = 0, int b = 0, ostream & out = cout) {
	int Times = 0, Check = 0;
	if (graph.empty()) {
		cout << "Graph is emty" << endl;
		return false;
	}
	int n = graph.size();
	if (a<1 || a>n || b<1 || b>n) {
		cout << "Nepravil'no vvedeno a ili b" << endl;
		return false;
	}
	a--; b--;

	bool result = graph[a][b];
	Graph prev = graph, cur(n, vector<int>(n, inf));
	for (int k = 2; k <= n - 1; k++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cur[i][j] = 0;
				for (int v = 0; v < n; v++) {
					Times++;
					int Sum = prev[i][v] == inf || graph[v][j] == inf ? inf : prev[i][v] + graph[v][j];
					if (Sum != inf && Sum > cur[i][j] || cur[i][j] == inf) {
						cur[i][j] = Sum;
					}
				}
				if (i == a && j == b && cur[i][j] != 0) {
					result = true;
				}
			}
		}
		if (k == count) {
			Check = Times;
		}
		if (k == count) {
			out << "Amount of edges: " << k << endl;
			out << cur << endl;
		}
		prev = cur;
	}
	cout << "���������� ��������(����������� ��������): " << Check << endl;
	//out << "����������� ���������� �������� �� ����� " << ++a << " � ����� " << ++b << ": ";
	return result;
}

bool Shimbell_new(const Graph& graph, int count, bool Print = true, int a = 0, int b = 0, ostream & out = cout) {
	int Times = 0, Check = 0;
	if (graph.empty()) {
		cout << "Graph is emty" << endl;
		return false;
	}
	int n = graph.size();
	if (a<1 || a>n || b<1 || b>n) {
		cout << "Nepravil'no vvedeno a ili b" << endl;
		return false;
	}
	a--; b--;

	bool result = graph[a][b];
	Graph prev = graph, cur(n, vector<int>(n));
	/*if (Print){
	out << "Amount of edges: 1" << endl;
	out << graph << endl;
	}*/
	for (int k = 2; k <= n - 1; k++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cur[i][j] = 0;
				for (int v = 0; v < n; v++) {
					Times++;
					int Sum = prev[i][v] == inf || graph[v][j] == inf ? inf : prev[i][v] + graph[v][j];
					if (Sum != inf && Sum < cur[i][j] || cur[i][j] == inf) {
						cur[i][j] = Sum;
					}
				}
				if (i == a && j == b && cur[i][j] != 0) {
					result = true;
				}
			}
		}
		if (k == count) {
			Check = Times;
		}
		if (k == count && Print) {
			out << "Amount of edges: " << k << endl;
			out << cur << endl;
		}
		prev = cur;
	}
	if (Print) {
		cout << "���������� ��������(����������� ��������): " << Check << endl;
	}
	//out << "����������� ���������� �������� �� ����� " << ++a << " � ����� " << ++b << ": ";
	return result;
}

vector<int> Dijkstra_new(const Graph& graph, int a, int b, bool print = true) {
	int Times = 0, Check = 0;
	int n = graph.size();
	if (a < 1 || a > n || b < 1 || b > n) {
		cout << "Graph vertex out of range." << endl;
		return vector<int>();
	}
	a--; b--;
	vector<int> Distances(n, inf), Prev(n, -1);
	vector<bool> Done(n);
	Distances[a] = 0;
	Prev[b] = a;
	while (true) {
		int V = -1;
		for (int i = 0; i < n; i++) {
			if (!Done[i] && (V == -1 || Distances[i] < Distances[V])) {
				V = i;
			}
		}
		if (V == b) {
			break;
		}
		Done[V] = true;
		for (int u = 0; u < n; u++) {
			if (Distances[V] + graph[V][u] < Distances[u]) {
				Distances[u] = Distances[V] + graph[V][u];
				Prev[u] = V;
			}
		}
	}
	//cout << graph.size() << endl;
	vector<int> Path;
	for (int i = b; i != -1; i = Prev[i]) {
		Path.push_back(i);
	}
	reverse(Path.begin(), Path.end());
	if (print) {
		cout << "Minimal path from " << a + 1 << " to " << b + 1 << ":\n";
		for (int i = 0; i < Path.size(); i++) {
			cout << Path[i] + 1 << " ";
		}
		cout << endl << "Length: " << Distances[b] << ".\n";
	}
	return Path;
}
//void Dijkstra(const Graph &graph, int a, int b){
//	clock_t t = clock();
//	clock_t finish;
//	//cout << t << endl;
//	a--;
//	b--;
//	int n = graph.size();
//	const int inf = 1000000000;
//	vector<int> Distances(n, inf);
//	vector<bool> Done(n);
//	Distances[a] = 0;
//	while (true){
//		int V = -1;
//		for (int i = 0; i < n; i++){
//			if (!Done[i] && (V == -1 || Distances[i] < Distances[V])){
//				V = i;
//			}
//		}
//		if (V == b){
//			break;
//		}
//		Done[V] = true;
//		for (int u = 0; u < n; u++){
//			if (graph[V][u] != 0 && Distances[V] + graph[V][u] < Distances[u]){
//				Distances[u] = Distances[V] + graph[V][u];
//			}
//		}
//
//	}
//	//t = clock();
//	int tmp = b;
//	tmp++;
//	cout << "���������� �� " << ++a << " �� " << tmp <<" ����� "<< Distances[b] << endl;
//	//finish = clock();
//	//cout << "Distance from a to b: " << Distances[b] << endl << "Time: " << (double)(finish - t) / CLOCKS_PER_SEC << endl;
//}
Graph To_undirected(const Graph& graph) {
	int n = graph.size();
	Graph result = graph;
	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			result[j][i] = graph[i][j];
		}
	}
	return result;
}

void Prim_algorithm(const Graph& graph) {
	Graph gr = To_undirected(graph);//������ ���� �����������������
	int Times = 0, Check = 0;
	int n = graph.size();
	Graph result(n, vector<int>(n, inf));
	vector<int> Distances(n, inf), Prev(n, -1);
	vector<bool> Done(n);
	int Sum = 0;
	int V = 0;
	Distances[V] = 0;
	for (int k = 0; k < n - 1; k++) {
		Done[V] = true;
		for (int u = 0; u < n; u++) {
			if (gr[V][u] < Distances[u]) {
				Distances[u] = gr[V][u];
				Prev[u] = V;
			}
		}
		V = -1;
		for (int i = 0; i < n; i++) {
			if (!Done[i] && (V == -1 || Distances[i] < Distances[V])) {
				V = i;
			}
		}
		result[V][Prev[V]] = result[Prev[V]][V] = gr[Prev[V]][V];//����������� ���������� �� V �� ����� �� ������ ������
		Sum += gr[Prev[V]][V];
	}
	cout << result << endl << "��������� ��� ������: " << Sum << endl << endl;

}

Graph Prim(const Graph& graph) {
	Graph gr = To_undirected(graph);
	int Times = 0, Check = 0;
	int n = graph.size();
	Graph result(n, vector<int>(n, inf));
	vector<int> Distances(n, inf), Prev(n, -1);
	vector<bool> Done(n);
	int Sum = 0;
	int V = 0;
	Distances[V] = 0;
	for (int k = 0; k < n - 1; k++) {
		Done[V] = true;
		for (int u = 0; u < n; u++) {
			if (gr[V][u] < Distances[u]) {
				Distances[u] = gr[V][u];
				Prev[u] = V;
			}
		}
		V = -1;
		for (int i = 0; i < n; i++) {
			if (!Done[i] && (V == -1 || Distances[i] < Distances[V])) {
				V = i;
			}
		}
		result[V][Prev[V]] = result[Prev[V]][V] = gr[Prev[V]][V];//����������� ���������� �� V �� ����� �� ������ ������
		Sum += gr[Prev[V]][V];
	}
	cout << result << endl << "��������� ��� ������: " << Sum << endl << endl;
	return result;
}

void Kruskal_algorithm(const Graph& graph) {
	int n = graph.size();
	Graph result(n, vector<int>(n, inf));
	int Sum = 0;
	vector<int> color(n);//�������������� � ����� �� ��������� ���-�� ������� ������ ���������
	vector<vector<int>> component(n, vector<int>(1));//���������� ���������(���-�� ������ ������� ������ � ����. ��������� ������� �������)
	//vector<int> sizeOfComponent(n,1);//������ ��������� ���������
	for (int i = 0; i < n; i++) {
		color[i] = i;
		component[i][0] = i;
	}
	Graph_edges edges = Matrix_to_Edges(graph);
	sort(edges.begin(), edges.end());
	for (int i = 0; i < edges.size(); i++) {
		int u = edges[i].second.first, v = edges[i].second.second;
		if (color[u] == color[v]) {
			continue;
		}
		result[u][v] = result[v][u] = edges[i].first;
		Sum += edges[i].first;
		u = color[u];
		v = color[v];
		if (component[u].size() < component[v].size()) {
			swap(u, v);
		}
		for (int k = 0; k < component[v].size(); k++) {
			color[component[v][k]] = u;//������������� � ���� ������� ����������
		}
		component[u].insert(component[u].end(), component[v].begin(), component[v].end());
	}
	cout << result << endl << "��������� ��� ������: " << Sum << endl << endl;
}

void Reachable(const Graph& graph, int a, int b, ostream& out = cout, bool Print = false) {
	a--;
	b--;
	bool result = graph[a][b];
	int n = graph.size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == a && j == b && graph[i][j] != 0) {
				result = true;
			}
		}
	}
	if (result) {
		cout << "���� ���� �� " << ++a << " � " << ++b << endl;
	}
	else { cout << "��� ���� �� " << ++a << " � " << ++b << endl; }
}

void Kirchhoff(Graph graph) {
	Graph graph1 = graph;
	for (int i = 0; i < graph.size(); i++) {
		for (int j = 0; j < graph.size(); j++) {
			if (graph1[i][j] == inf) { graph1[i][j] = 0; }
			if (graph1[i][j] != 0) { graph1[i][j] = 1; }
		}
	}
	Graph Kirchhoff(graph.size());
	Graph Alg(graph.size() - 1, vector<int>(graph.size() - 1));
	Kirchhoff = graph1;
	vector <int> deg;
	for (int i = 0; i < graph.size(); i++) {
		for (int j = 0; j < i; j++) {
			Kirchhoff[i][j] = Kirchhoff[j][i];
		}
	}
	for (int i = 0; i < graph.size(); i++) {
		deg.push_back(0);
		for (int j = 0; j < graph.size(); j++) {
			deg[i] += Kirchhoff[i][j];
		}
	}
	for (int i = 0; i < graph.size(); i++) {
		for (int j = 0; j < graph.size(); j++) {
			if (i == j) Kirchhoff[i][j] = deg[i];
			else if (Kirchhoff[i][j]) Kirchhoff[i][j] = -1;
		}
	}
	//cout << "Kirchhoff's matrix: \n";
	//Kirchhoff.printGraph();
	cout << endl;
	for (int i = 0; i < Alg.size(); i++) {
		for (int j = 0; j < Alg.size(); j++) {
			Alg[i][j] = Kirchhoff[i + 1][j + 1];
		}
	}
	cout << "the number of spanning trees is " << determinant(Alg) << endl;
}

vector<int> PruferCode(Graph graph) {
	Graph tree(graph.size(), vector<int>(graph.size()));
	Graph graph1 = Prim(graph);
	vector <int> deg = ddeg(graph1);//������� ������ �����
	vector <int> ver;
	vector <int> prufer;
	vector <int> decode;
	bool flag = true;
	bool flag2 = true;
	int i = 0;
	int j = 0;
	//	To_undirected(graph);

	for (int quu = 0; quu < graph.size(); quu++) {
		ver.push_back(quu);//������ ���� ������
	}
	while (ver.size() > 2) {//���������� ���� �������
		while (deg[i] != 1) i++;//���� ����
		//	while (!graph1[ver[i]][ver[j]]) j++;
		for (int q = 0; q < ver.size(); q++) {
			if (i != j && graph1[ver[i]][ver[j]] != inf && graph1[ver[i]][ver[j]] != 0) {//���� ������� ������� � ���� ������
				break;
			}
			j++;
		}
		prufer.push_back(ver[j]);//��������� � ��� ������� ����� ������� � ������ �������
		/*for (int k = 0; k < graph.size(); k++){
		if (graph[k][i] != inf && graph[k][i] != 0){
		deg[k]--;
		}
		}*/
		deg[j]--;//�������� ������� ����������� ������� �� 1
		deg.erase(deg.begin() + i);//������� �������� ������� ����� �� ������� �������� ������ �����
		ver.erase(ver.begin() + i);//������� ��������� ���� �� ������� ���� ������ �����
		i = 0; j = 0;
	}
	cout << "��� �������: ";
	for (int it = 0; it < prufer.size(); it++) {
		cout << prufer[it] + 1;
	}
	//�������������
	cout << endl;
	ver.clear();
	for (int quu = 0; quu < graph.size(); quu++) {
		ver.push_back(quu);//��������� ������ ���� ������ �����
	}
	for (int k = 0; k < ver.size(); k++) {//���� ������� � ����������� �������, �� ������������ � ���� �������
		for (int q = 0; q < prufer.size(); q++) {
			if (k == prufer[q]) {
				flag = false;
				break;
			}
		}
		if (flag) decode.push_back(k);//���������� ��������������� ��� ������� �������������� � ���� �������
		flag = true;
	}
	while (prufer.size() != 0) {
		tree[decode[0]][prufer[0]] = 1;//������� ����� ����� ����������� �������� �� ���� ������� � ����������� �������� ����� �� ������������ � ���� �������
		tree[prufer[0]][decode[0]] = 1;//� � �������� ������� ����
		i = decode[0];
		j = prufer[0];
		prufer.erase(prufer.begin());//������� ������ ������� �� �������
		decode.erase(decode.begin());//������� ������ ������� �� ������� ������ �� �������� � ��� �������
		for (int q = 0; q < prufer.size(); q++) {
			if (j == prufer[q]) {
				flag = false;
				break;
			}
		}
		if (flag) {
			for (int qu = 0; qu < decode.size(); qu++) {
				if (j < decode[qu]) {
					decode.insert(decode.begin() + qu, j);
					flag2 = false;
					break;
				}
			}
			if (flag2) decode.push_back(j);
			flag2 = true;
		}
		flag = true;
	}
	tree[j][ver.size() - 1] = 1;//��������� ��������� 2 ���������� ������� ������
	tree[ver.size() - 1][j] = 1;//
	cout << "�������������� �����: \n";
	cout << tree << endl;
	return prufer;
}

void InfToZeros(Graph& graph)
{
	int n = graph.size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (graph[i][j] == inf) {
				graph[i][j] = 0;
			}
			else if (graph[i][j] == 0) {
				graph[i][j] = 1;
			}
		}
	}
}

void ZerosToInf(Graph& graph)
{
	int n = graph.size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (graph[i][j] == 0) {
				graph[i][j] = inf;
			}
		}
	}
}
bool isEilerGraph(const vector<vector<int>>& graph) {//��������� ���� �� �����������, ���� �� ������� � ��������� ���������?
	int chetn = 0;
	Graph graph1 = graph;
	for (int i = 0; i < graph.size(); i++) {
		for (int j = 0; j < graph.size(); j++) {
			if (graph1[i][j] == inf) { graph1[i][j] = 0; }
			if (graph1[i][j] != 0) { graph1[i][j] = 1; }
		}
	}
	for (int i = 0; i < graph.size(); i++)
	{
		for (int j = 0; j < graph.size(); j++)
		{
			if (graph1[i][j])
				chetn++;
		}
		if (chetn & 1)
			return false;
		chetn = 0;
	}
	return true;
}

bool isNull(const Graph& graph) {//���� �� ��� ����� � �������?
	int n = graph.size();
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (graph[i][j] != inf && graph[i][j] != 0) return false;
	return true;
}

bool Svyazn(const Graph& matr, int tv, int in) {//�������� ������-�������� � ������������ �����������. ����������������� ���� ������, ���� �� ������ ������� ����� ������� �� ��� ���������.
	vector<vector<int>> temp(matr);//� �������� ������ ��� ��� ���� ���� ���� �� ���� �� ����
	int n = matr.size();
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				if (i != k && i != j && temp[j][i] != 0 && temp[i][k] != 0 && (min(temp[j][k], (temp[j][i] + temp[i][k])) != temp[j][k]))//
					temp[j][k] = temp[j][i] + temp[i][k];
	for (int i = 1; i < n; i++)
		if (!temp[0][i] && noway(matr, i, tv, in)) return false;//���� �� �����-�� ������� (����� ���� �����) ���� ���, ���� ��������
	return true;
}

vector<int> Eiler(vector<vector<int>> matr) {//������� ����� �� ���� ������, ��� ��� ����� �� �� ��������, � �� ������
	int n = matr.size();
	int iter = 0;
	vector<int> rezult;
	int tv = 0, tmp;
	bool b = true;
	rezult.push_back(1);
	while (!isNull(matr)) {
		b = false;
		for (int i = n - 1; i >= 0; i--) {
			iter++;
			if (matr[tv][i]) {
				tmp = matr[tv][i];
				matr[i][tv] = matr[tv][i] = 0;
				if (!b && i != 0 && Svyazn(matr, tv, i) || b) {//
					rezult.push_back(i + 1);
					tv = i;
					break;
				}//
				else
					matr[i][tv] = matr[tv][i] = tmp;//
			}
			if (i == 0) {//
				i = n;//
				b = true;
			}
		}
	}
	return rezult;
}

//bool gham(Graph g, int curr, int n, vector<int> &Path, vector<bool> &Visited)
//{
//	Path.push_back(curr);
//	if (Path.size() == n)
//	if (g[Path[0]][Path.back()])
//		return true;
//	else
//	{
//		Path.pop_back();
//		return false;
//	}
//	Visited[curr] = true;
//
//	//min
//	vector<int> minp;
//	vector<int> stock;
//	int count = 0;
//	for (int i = 0; i < n; i++)
//	{
//		if (g[curr][i]>0)
//		{
//			minp.push_back(g[curr][i]);
//			stock.push_back(i);
//			count++;
//		}
//	}
//	for (int i = 0; i < count; i++)
//	{
//		for (int j = 0; j < count - 1; j++)
//		{
//			if (minp[j]>minp[j + 1])
//			{
//				swap(minp[j], minp[j + 1]);
//				swap(stock[j], stock[j + 1]);
//			}
//		}
//	}
//
//	//end min
//
//	for (int nxt = 0; nxt < count; ++nxt)
//	{
//
//		if (g[curr][stock[nxt]] >0 && !Visited[stock[nxt]])
//		if (gham(g, stock[nxt], n, Path, Visited))
//			return true;
//	}
//	Visited[curr] = false;
//	Path.pop_back();
//
//	return false;
//}
//
//void hamilton(Graph graph){
//	Graph Hamilton = graph;
//	int tmp = Hamilton.size();//Hamilton.graph.size()
//	Graph tem = Hamilton;
//
//	Graph graph1 = graph;
//	for (int i = 0; i < graph.size(); i++){
//		for (int j = 0; j < graph.size(); j++) {
//			if (graph1[i][j] == inf){ graph1[i][j] = 0; }
//			if (graph1[i][j] != 0){ graph1[i][j] = 1; }
//		}
//	}
//	vector<int> pjj;
//	vector < vector < int >> arrways;
//	vector <bool> visited(Hamilton.size(), false);
//
//
//
//	for (int i = 0; i < graph.size(); i++)
//	{
//		if (gham(tem, i, Hamilton.size(), pjj, visited))
//		{
//
//			arrways.push_back(pjj);
//			pjj.clear();
//			for (int j = 0; j < visited.size(); j++)
//				visited[j] = false;
//		}
//	}
//
//	///end second
//	vector <bool> point;
//	vector <int> ver;
//	int temp;
//	int temp1 = graph.size() - 1;
//	bool flag;
//
//
//	cout << Hamilton << endl;
//	
//	//ver = pjj;
//
//
//	//int c1 = 0;
//	cout << endl;
//	cout << tem << endl;
//	//bool a = hamCycle(tem);
//	//output hamiltonian cycles
//	vector<vector<int>> cycles;
//	vector<int> costs;
//	for (int j = 0; j < arrways.size(); j++)
//	{
//		//cout << "Hamilton's graph is: ";
//		int c1 = 0;
//		for (size_t i = 1; i < arrways[j].size(); i++){
//			//cout << arrways[j][i] << " ";
//			c1 += tem[arrways[j][i - 1]][arrways[j][i]];
//		}
//		int last = arrways[j][arrways[j].size() - 1];
//		c1 += tem[last][arrways[j][0]];
//		bool rep = true;
//
//		cycles.push_back(arrways[j]);
//		costs.push_back(c1);
//
//		//cout << arrways[j][0] << " Cost: " << c1 << endl;
//	}
//	cout << endl << "Hamilton's cycles:" << endl;
//	for (int c = 0; c < costs.size(); c++)
//	{
//		for (size_t i = 0; i < cycles[c].size(); i++){
//			cout << cycles[c][i] << " ";
//
//		}
//		cout << cycles[c][0] << " Cost: " << costs[c] << endl;;
//	}
//	int minc = costs[0];
//	for (int i = 1; i < costs.size(); i++) if (costs[i] < minc)minc = costs[i];
//	cout << "Minimal Cost: " << minc << endl;
//
//}

int Salesman(const Graph& graph, vector<int>& cycle, int v = 0)//����� ����������� ����� ������������ �����
{
	static vector<int> path;
	static vector<bool> busy;

	int n = graph.size();
	busy.resize(n);
	int result = inf;

	path.push_back(v);//�������� � 0-�� �������
	if (path.size() == n) {
		cycle = path;
		int w = graph[v][path[0]];
		if (w < inf) {
			int Length = 0;
			for (int i = 0; i < n; i++) {
				cout << cycle[i] + 1 << "->";
				Length += graph[cycle[i]][cycle[(i + 1) % n]];
			}
			cout << cycle[0] + 1 << "  ";
			cout << "���: " << Length << endl;
		}
		result = w;
	}
	else {
		busy[v] = true;
		for (int i = 0; i < n; ++i) {
			if (!busy[i] && graph[v][i] < inf) {//���� ����� ������������ �������, ��������� ������� ���������� �� ���
				vector<int> cycle0;
				int length = Salesman(graph, cycle0, i);
				if (length < inf && graph[v][i] + length < result) {
					result = graph[v][i] + length;
					cycle = cycle0;
				}
			}
		}
		busy[v] = false;
	}
	path.pop_back();

	return result;
}

Graph Hamiltonian(const Graph& graph)//���������� ����� � ������������
{
	Graph graph1 = graph;
	int n = graph.size();
	vector<bool> done(n);
	int v = 0;
	for (int count = 1; count < n; ++count) {
		done[v] = true;
		int i;
		for (i = 0; i < n && (done[i] || graph1[v][i] == inf); ++i);
		if (i == n) {
			for (i = 0; done[i]; ++i);
			graph1[v][i] = graph1[i][v] = 1;
		}
		v = i;
	}
	if (graph1[v][0] == inf) {
		graph1[v][0] = graph1[0][v] = 1;
	}
	return graph1;
}

void Hamilton(const Graph& graph)
{
	vector<int> cycle;
	int length = Salesman(graph, cycle);
	if (length == inf) {
		cout << "\n���� �� �����������. �������� � ������������:\n\n";
		Graph ham = Hamiltonian(graph);
		cout << ham;
		length = Salesman(ham, cycle);
	}
	cout << "\n���������� ����������� ����:\n";
	for (int i = 0; i < cycle.size(); ++i) {
		cout << cycle[i] + 1 << " -> ";
	}
	cout << 1 << "\n����� �����: " << length << endl << endl;
}

Network GenerateNetwork(const Graph& graph) {//��������� �������. ������. � ����������
	int n = graph.size();
	int Sum = 0;//���-�� �����
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (graph[i][j] != inf) {
				Sum++;
			}
		}
	}
	Network network(n, vector<e>(n));
	vector<int> Capacity = GenerateBinomial(Sum, 99, 0.05), Price = GenerateBinomial(Sum, 99, 0.2);
	int index = 0;
	for (int i = 0; i < graph.size(); i++) {
		for (int j = 0; j < graph.size(); j++) {
			if (graph[i][j] != inf) {
				network[i][j].c = Capacity[index];
				network[i][j].a = Price[index];
				index++;
			}
		}
	}
	return network;
}

int Fulkerson(const Network& network, int s, int t) {
	int n = network.size();
	int f = 0;
	Graph ResidualCapacity(n, vector<int>(n));//������� ����� � ���������� ���������� ������������
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			ResidualCapacity[i][j] = network[i][j].c;
		}
	}
	Graph forDijkstra = ResidualCapacity;
	while (true) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				forDijkstra[i][j] = ResidualCapacity[i][j] > 0 ? 1 : inf;//���� ���� �������. �������. ������. ������ ����� ������������ ����� [i][j]  ��� ������ ����. ����
			}
		}
		vector<int> path = Dijkstra_new(forDijkstra, s, t, false);//���� ��������� ����
		if (forDijkstra[path[0]][path[1]] == inf) {
			break;
		}
		for (int i = 0; i < path.size() - 1; i++) {
			int u = path[i], v = path[i + 1];
			ResidualCapacity[u][v]--;//��������� �������. �������. ������.
			ResidualCapacity[v][u]++;//��������� ����������� ��������� ����� �� ����� [u][v]
		}
		f++;
	}
	Graph FlowMatrix(n, vector<int>(n));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			FlowMatrix[i][j] = network[i][j].c - ResidualCapacity[i][j];
		}
	}
	cout << "������� ������������� ����� ��� ������������ ������: " << endl << FlowMatrix << endl << endl;
	return f;

}

int FindMinCostFlow(const Network& network, int s, int t, int flow) {
	int n = network.size(), cost = 0;
	Graph FlowMatrix(n, vector<int>(n));
	Graph ResidualCapacity(n, vector<int>(n));//������� ����� � ���������� ���������� ������������
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			ResidualCapacity[i][j] = network[i][j].c;
		}
	}
	Graph MatrixOfPrices = ResidualCapacity;
	for (int f = 0; f < flow; f++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				MatrixOfPrices[i][j] = inf;
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (ResidualCapacity[i][j] > 0) {//���� ���� ��� ��������. ������� �����������
					MatrixOfPrices[i][j] = network[i][j].c > 0 ? network[i][j].a : -network[i][j].a;//���� ���� �������. �������. ������. ������ ����� ������������ ����� [i][j]  ��� ������ ����. ����, 0 � ������� ���������� ����� ���, ���� ��� ������ ������� �����
				}
			}
		}
		vector<int> path = Bellman_Ford(MatrixOfPrices, s, t, false);//���� ��������� ����
		for (int i = 0; i < path.size() - 1; i++) {
			int u = path[i], v = path[i + 1];
			ResidualCapacity[u][v]--;//��������� �������. �������. ������.
			ResidualCapacity[v][u]++;//��������� ����������� ��������� ����� �� ����� [u][v]
			cost += MatrixOfPrices[u][v];

		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				FlowMatrix[i][j] = network[i][j].c - ResidualCapacity[i][j];
			}
		}
		cout << FlowMatrix << endl << endl;
	}

	/*for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			FlowMatrix[i][j] = network[i][j].c - ResidualCapacity[i][j];
		}
	}*/
	cout << "������� ������������� ����� ��� 2/3 �� ������������� ������: " << endl << FlowMatrix << endl << endl;
	//cout << "������� ������������� �����:" << endl;

	return cost;
}
}