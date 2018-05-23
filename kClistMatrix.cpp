/*
Info:
Feel free to use these lines as you wish.
This program iterates over all k-cliques. It is highly scallable.

To compile:
"gcc kClistDens.c -O9 -o kClistDens -fopenmp".

To execute:
"./kClistDens p k edgelist.txt".
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
k is the size of the k-cliques
p is the number of threads
Will print the number of k-cliques and the density of the found kclique densest.
*/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <unordered_map>
#include <vector>


#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

class Clique_Matrix {
	public:
		Clique_Matrix(unsigned);
		unsigned vector_length;
		std::vector<std::unordered_map<unsigned, int64_t>> clique_mat;
		void add_edge(unsigned, unsigned);
};

Clique_Matrix::Clique_Matrix(unsigned number_of_nodes) {
	vector_length = number_of_nodes;
	clique_mat = std::vector<std::unordered_map<unsigned, int64_t>>(number_of_nodes);
};

void Clique_Matrix::add_edge(unsigned p, unsigned q) {
	clique_mat[p][q]++;
}

typedef struct {
	unsigned s;
	unsigned t;
} edge;

typedef struct {
	unsigned node;
	unsigned deg;
} nodedeg ;

typedef struct {
	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge *edges;//list of edges
	unsigned *rank;//ranking of the nodes according to degeneracy ordering
	//unsigned *map;//oldID newID correspondance NOT USED IN THIS VERSION
	unsigned *node_map; //map from new id to old id
} edgelist;

typedef struct {
	unsigned n;
	unsigned e;
	edge *edges;//ading this again here: TO IMPROVE
	unsigned *cd;//cumulative degree: (starts with 0) length=n+1
	unsigned *adj;//truncated list of neighbors
	unsigned core;//core value of the graph
} graph;

typedef struct {
	unsigned *n;//n[l]: number of nodes in G_l
	unsigned **d;//d[l]: degrees of G_l
	unsigned *adj;//truncated list of neighbors
	unsigned char *lab;//lab[i] label of node i
	unsigned **nodes;//sub[l]: nodes in G_l
	unsigned core;
} subgraph;

void free_graph(graph *g){
	free(g->cd);
	free(g->adj);
	free(g);
}

void free_subgraph(subgraph *sg, unsigned char k){
	unsigned char i;
	free(sg->n);
	for (i=1;i<k;i++){
		free(sg->d[i]);
		free(sg->nodes[i]);
	}
	free(sg->d);
	free(sg->nodes);
	free(sg->lab);
	free(sg->adj);
	free(sg);
}


//Compute the maximum of three unsigned integers.
inline unsigned int max3(unsigned int a,unsigned int b,unsigned int c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

edgelist* readedgelist(char* input){
	unsigned e1=NLINKS;
	edgelist *el=(edgelist*)malloc(sizeof(edgelist));
	FILE *file;

	el->n=0;
	el->e=0;
	file=fopen(input,"r");
	el->edges=(edge*)malloc(e1*sizeof(edge));
	unsigned w = 1;
	unsigned s = 1;
	unsigned t = 1;

        fscanf(file, "%u %u %u", %s, %t, %w);
	while (fscanf(file,"%u %u %u", &s, &t, &w)==3) {//Add one edge
		if (s < t) {
			(el->edges[el->e].s) = s;
			(el->edges[el->e].t) = t;
			el->n=max3(el->n,el->edges[el->e].s,el->edges[el->e].t);
			el->e++;
			if (el->e==e1) {
				e1+=NLINKS;
				el->edges=(edge*)realloc(el->edges,e1*sizeof(edge));
			}
		}
	}
	fclose(file);
	el->n++;

	el->edges=(edge*)realloc(el->edges,el->e*sizeof(edge));
	el->node_map=(unsigned*)malloc((el->n)*sizeof(unsigned));


	return el;
}

void relabel(edgelist *el){
	// for (int m=0; m<el->n; m++) {
	// 	printf("node: %u  rank: %u\n", m, el->rank[m]);
	// }
	unsigned i, source, target, tmp;
	el->n=0;
		//FILE* file=fopen("debug.txt","w");
	for (i=0;i<el->e;i++) {
		source=el->rank[el->edges[i].s];
		target=el->rank[el->edges[i].t];

		// printf("old src: %u  new src: %u\n", el->edges[i].s, source );
		// printf("old t: %u  new t: %u\n", el->edges[i].t, target );
		el->node_map[source] =  el->edges[i].s;
		el->node_map[target] = el->edges[i].t;

		if (source<target){
			tmp=source;
			source=target;
			target=tmp;
		}
		if (source+1>el->n){
			el->n=source+1;
		}


		el->edges[i].s=source;
		el->edges[i].t=target;

		//fprintf(file,"%u %u\n",source,target);
	}
	//fclose(file);

}

///// CORE ordering /////////////////////

typedef struct {
	unsigned key;
	unsigned value;
} keyvalue;

typedef struct {
	unsigned n_max;	// max number of nodes.
	unsigned n;	// number of nodes.
	unsigned *pt;	// pointers to nodes.
	keyvalue *kv; // nodes.
} bheap;


bheap *construct(unsigned n_max){
	unsigned i;
	bheap *heap=(bheap*)malloc(sizeof(bheap));

	heap->n_max=n_max;
	heap->n=0;
	heap->pt=(unsigned *)malloc(n_max*sizeof(unsigned));
	for (i=0;i<n_max;i++) heap->pt[i]=-1;
	heap->kv=(keyvalue*)malloc(n_max*sizeof(keyvalue));
	return heap;
}

void swap(bheap *heap,unsigned i, unsigned j) {
	keyvalue kv_tmp=heap->kv[i];
	unsigned pt_tmp=heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key]=heap->pt[heap->kv[j].key];
	heap->kv[i]=heap->kv[j];
	heap->pt[heap->kv[j].key]=pt_tmp;
	heap->kv[j]=kv_tmp;
}

void bubble_up(bheap *heap,unsigned i) {
	unsigned j=(i-1)/2;
	while (i>0) {
		if (heap->kv[j].value>heap->kv[i].value) {
			swap(heap,i,j);
			i=j;
			j=(i-1)/2;
		}
		else break;
	}
}

void bubble_down(bheap *heap) {
	unsigned i=0,j1=1,j2=2,j;
	while (j1<heap->n) {
		j=( (j2<heap->n) && (heap->kv[j2].value<heap->kv[j1].value) ) ? j2 : j1 ;
		if (heap->kv[j].value < heap->kv[i].value) {
			swap(heap,i,j);
			i=j;
			j1=2*i+1;
			j2=j1+1;
			continue;
		}
		break;
	}
}

void insert(bheap *heap,keyvalue kv){
	heap->pt[kv.key]=(heap->n)++;
	heap->kv[heap->n-1]=kv;
	bubble_up(heap,heap->n-1);
}

void update(bheap *heap,unsigned key){
	unsigned i=heap->pt[key];
	if (i!=-1){
		((heap->kv[i]).value)--;
		bubble_up(heap,i);
	}
}

keyvalue popmin(bheap *heap){
	keyvalue min=heap->kv[0];
	heap->pt[min.key]=-1;
	heap->kv[0]=heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key]=0;
	bubble_down(heap);
	return min;
}

//Building the heap structure with (key,value)=(node,degree) for each node
bheap* mkheap(unsigned n,unsigned *v){
	unsigned i;
	keyvalue kv;
	bheap* heap=construct(n);
	for (i=0;i<n;i++){
		kv.key=i;
		kv.value=v[i];
		insert(heap,kv);
	}
	return heap;
}

void freeheap(bheap *heap){
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

//computing degeneracy ordering and core value
void ord_core(edgelist* el){
	unsigned i,j,r=0,n=el->n,e=el->e;
	keyvalue kv;
	bheap *heap;

	unsigned *d0=(unsigned*)calloc(el->n,sizeof(unsigned));
	unsigned *cd0=(unsigned *)malloc((el->n+1)*sizeof(unsigned));
	unsigned *adj0=(unsigned *)malloc(2*el->e*sizeof(unsigned));
	for (i=0;i<e;i++) {
		d0[el->edges[i].s]++;
		d0[el->edges[i].t]++;
	}
	cd0[0]=0;
	for (i=1;i<n+1;i++) {
		cd0[i]=cd0[i-1]+d0[i-1];
		d0[i-1]=0;
	}
	for (i=0;i<e;i++) {
		adj0[ cd0[el->edges[i].s] + d0[ el->edges[i].s ]++ ]=el->edges[i].t;
		adj0[ cd0[el->edges[i].t] + d0[ el->edges[i].t ]++ ]=el->edges[i].s;
	}

	heap=mkheap(n,d0);

	el->rank=(unsigned *)malloc(n*sizeof(unsigned));
	for (i=0;i<n;i++){
		kv=popmin(heap);
		el->rank[kv.key]=n-(++r);
		for (j=cd0[kv.key];j<cd0[kv.key+1];j++){
			update(heap,adj0[j]);
		}
	}
	freeheap(heap);
	free(d0);
	free(cd0);
	free(adj0);
}

//////////////////////////
//Building the special graph
graph* mkgraph(edgelist *el){
	unsigned i,max;
	unsigned *d;
	graph* g=(graph*)malloc(sizeof(graph));

	d=(unsigned *)calloc(el->n,sizeof(unsigned));

	for (i=0;i<el->e;i++) {
		d[el->edges[i].s]++;
	}

	g->cd=(unsigned *)malloc((el->n+1)*sizeof(unsigned));
	g->cd[0]=0;
	max=0;
	for (i=1;i<el->n+1;i++) {
		g->cd[i]=g->cd[i-1]+d[i-1];
		max=(max>d[i-1])?max:d[i-1];
		d[i-1]=0;
	}
	//printf("core value (max truncated degree) = %u\n",max);

	g->adj=(unsigned *)malloc(el->e*sizeof(unsigned));

	for (i=0;i<el->e;i++) {
		g->adj[ g->cd[el->edges[i].s] + d[ el->edges[i].s ]++ ]=el->edges[i].t;
	}

	free(d);
	g->core=max;
	g->n=el->n;

	free(el->rank);
	g->edges=el->edges;
	g->e=el->e;
	//free(el);
	//		printf("el2=%u\n",el->e);
	return g;
}


subgraph* allocsub(graph *g,unsigned char k){
	unsigned i;
	subgraph* sg=(subgraph *)malloc(sizeof(subgraph));
	sg->n=(unsigned *)calloc(k,sizeof(unsigned));
	sg->d=(unsigned **)malloc(k*sizeof(unsigned*));
	sg->nodes=(unsigned **)malloc(k*sizeof(unsigned*));
	for (i=1;i<k;i++){/////////
		sg->d[i]=(unsigned *)malloc(g->core*sizeof(unsigned));
		sg->nodes[i]=(unsigned *)malloc(g->core*sizeof(unsigned));
	}
	sg->lab=(unsigned char *)calloc(g->core,sizeof(unsigned char));
	sg->adj=(unsigned *)malloc(g->core*g->core*sizeof(unsigned));
	sg->core=g->core;
	return sg;
}


unsigned *old=NULL,*newN=NULL;//to improve
#pragma omp threadprivate(newN,old)

void mksub(graph* g,edge ed,subgraph* sg,unsigned char k){
	unsigned i,j,l,x,y;
	unsigned u=ed.s,v=ed.t;

	if (old==NULL){
		newN=(unsigned *)malloc(g->n*sizeof(unsigned));
		old=(unsigned *)malloc(g->core*sizeof(unsigned));
		for (i=0;i<g->n;i++){
			newN[i]=-1;
		}
	}

	for (i=0;i<sg->n[k-1];i++){
		sg->lab[i]=0;
	}

	for (i=g->cd[v];i<g->cd[v+1];i++){
		newN[g->adj[i]]=-2;
	}

	j=0;
	for (i=g->cd[u];i<g->cd[u+1];i++){
		x=g->adj[i];
		if (newN[x]==-2){
			newN[x]=j;
			old[j]=x;
			sg->lab[j]=k-2;
			sg->nodes[k-2][j]=j;
			sg->d[k-2][j]=0;//new degrees
			j++;
		}
	}

	sg->n[k-2]=j;

	for (i=0;i<sg->n[k-2];i++){//reodering adjacency list and computing new degrees
		x=old[i];
		for (l=g->cd[x];l<g->cd[x+1];l++){
			y=g->adj[l];
			j=newN[y];
			if (j<-2){
				sg->adj[sg->core*i+sg->d[k-2][i]++]=j;
			}
		}
	}

	for (i=g->cd[v];i<g->cd[v+1];i++){
		newN[g->adj[i]]=-1;
	}
}


unsigned long long *ckdeg_p,*ckdeg;
unsigned * ck_buf;
unsigned *ck_p;
Clique_Matrix * ck_m;
Clique_Matrix *clique_matrix;
int pos;
#pragma omp threadprivate(ck_m, ckdeg_p,ck_p,ck_buf,pos)

	//after the threads are done, combine individual clique matrices to overall one
	void combine_to_global_matrix() {
		#pragma omp critical (mat)
		for (int p = 0; p < ck_m->vector_length; p++) {//vector iterator
			for (std::pair<unsigned, int64_t> element : ck_m->clique_mat[p]) {//map iterator **MAKE THIS FASTER BY ONLY ADDING ONES THERE WITH ITER
					clique_matrix->clique_mat[p][element.first]+= element.second;
			}
		}
	}

	void add_to_thread_matrix(unsigned k) {
		for(int i=0; i<pos; i=i+k) {  //<pos
			for(int p=0; p < k; p++) {
				for(int q=p+1;q<k;q++) {
					ck_m->add_edge(ck_buf[i+p], ck_buf[i+q]);
					ck_m->add_edge(ck_buf[i+q], ck_buf[i+p]);
				}
			}
		}
	}

	void add_clique_to_buf(unsigned k, unsigned * clique) {
		for (int i = 0; i < k; i++) {
			ck_buf[pos] = clique[i];
			pos++;
		}
		if (pos >= 10000*k) {
			add_to_thread_matrix(k);
			pos=0;
		}
	}

void allocglobal(graph *g,unsigned k, unsigned number_of_nodes){
        #pragma omp parallel
        {
                ck_p=(unsigned *)calloc(k,sizeof(unsigned));
                ckdeg_p=(unsigned long long *)calloc(g->n,sizeof(unsigned long long));
								ck_buf=(unsigned *)calloc(10000*k, sizeof(unsigned));
								pos = 0;
								ck_m = new Clique_Matrix(number_of_nodes);
        }
        ckdeg=(unsigned long long *)calloc(g->n,sizeof(unsigned long long));
				clique_matrix = new Clique_Matrix(number_of_nodes);
}

void kclique_thread(unsigned char kmax, unsigned char l, subgraph *sg, unsigned long long *n, unsigned * node_map) {
	unsigned i,j,k,end,u,v,w;

	if (kmax==3){//can be improved
		for(i=0; i<sg->n[1]; i++){//list all nodes
			ckdeg_p[old[sg->nodes[1][i]]]++;
			ckdeg_p[ck_p[1]]++;
			ckdeg_p[ck_p[2]]++;
			(*n)++;//listing here!!!
			//ADD EDGELIST AND PRINT MAPPED NODES
			unsigned kclique[3] = {node_map[ck_p[1]], node_map[ck_p[2]], node_map[old[sg->nodes[1][i]]]};
			add_clique_to_buf(3, kclique);
			//printf("ADD CLIQUE %u,%u,%u\n", node_map[ck_p[1]], node_map[ck_p[2]], node_map[old[sg->nodes[1][i]]]);

		}
		for (i=0;i < (sizeof (ckdeg_p) /sizeof (ckdeg_p[0]));i++) {
		}
		return;
	}

	if(l==2){

		for(i=0; i<sg->n[2]; i++){//list all edges
			u=sg->nodes[2][i];
			end=u*sg->core+sg->d[2][u];
			for (j=u*sg->core;j<end;j++) {
				unsigned * kclique = (unsigned *)malloc(kmax*sizeof(unsigned));
				int count = 0;
				//directly below is edge u,v
				//insert u, v and v,  u
				kclique[count] = node_map[old[u]];
				count++;
				kclique[count]= node_map[old[sg->adj[j]]];
				count++;
				ckdeg_p[old[sg->adj[j]]]++;
				ckdeg_p[old[u]]++;
				for (l=2;l<kmax;l++){//ok to use l here :)
					ckdeg_p[ck_p[l]]++;
					kclique[count] = node_map[ck_p[l]];
					count++;
				}
				(*n)++;//listing here!!!
				add_clique_to_buf(kmax, kclique);
				for (int m = 0; m <kmax; m++){
					//printf("%u, ",kclique[m]);
				}
				//printf("\n");
			}

		}

		return;
	}

	for(i=0; i<sg->n[l]; i++){
		u=sg->nodes[l][i];
		ck_p[l-1]=old[u];
		//printf("%u %u\n",i,u);
		sg->n[l-1]=0;
		end=u*sg->core+sg->d[l][u];
		for (j=u*sg->core;j<end;j++){//relabeling nodes and forming U'.
			v=sg->adj[j];
			if (sg->lab[v]==l){
				sg->lab[v]=l-1;
				sg->nodes[l-1][sg->n[l-1]++]=v;
				sg->d[l-1][v]=0;//new degrees
			}
		}
		for (j=0;j<sg->n[l-1];j++){//reodering adjacency list and computing new degrees
			v=sg->nodes[l-1][j];
			end=sg->core*v+sg->d[l][v];
			for (k=sg->core*v;k<end;k++){
				w=sg->adj[k];
				if (sg->lab[w]==l-1){
					sg->d[l-1][v]++;
				}
				else{
					sg->adj[k--]=sg->adj[--end];
					sg->adj[end]=w;
				}
			}
		}

		kclique_thread(kmax,l-1, sg, n, node_map);

		for (j=0;j<sg->n[l-1];j++){//restoring labels
			v=sg->nodes[l-1][j];
			sg->lab[v]=l;
		}

	}
}

unsigned long long kclique_main(unsigned char k, graph *g, unsigned * node_map) {
	unsigned i;
	unsigned long long n=0;
	subgraph *sg;
	#pragma omp parallel private(sg,i) reduction(+:n)
	{
		sg=allocsub(g,k);

		#pragma omp for schedule(dynamic, 1) nowait
		for(i=0; i<g->e; i++){
			ck_p[k-1]=g->edges[i].s;
			ck_p[k-2]=g->edges[i].t;
			mksub(g,g->edges[i],sg,k);
			kclique_thread(k,k-2, sg, &n,node_map);

		}
		add_to_thread_matrix(k);
		pos=0;
		combine_to_global_matrix();
		free_subgraph(sg,k);
		#pragma omp single
		{
		bzero(ckdeg,g->n*sizeof(unsigned long long));
		}

		#pragma omp barrier //is it necessary???

		#pragma omp critical
		{
			for(i=0; i<g->n; i++){
				ckdeg[i]+=ckdeg_p[i];
			}
			bzero(ckdeg_p,g->n*sizeof(unsigned long long));
		}

	}
	return n;
}

void rmnodes(bool *rm,edgelist* el){
	unsigned long long i,r=0;
	FILE* file=fopen("debug.txt","w");
		for (i=0;i<el->e;i++){
				if (((rm[el->edges[i].s]==1) || (rm[el->edges[i].t]==1)) == 0){
					r++;
									fprintf(file,"%u %u\n",el->edges[i].s,el->edges[i].t);
				}
		}
	//printf("el0=%llu\n",r);
	for (i=0;i<el->e;i++){
//printf("%llu\n",i);
		if ((rm[el->edges[i].s]==1) || (rm[el->edges[i].t]==1)){
			el->edges[i--]=el->edges[--(el->e)];
		}
	}
			//printf("el1=%u\n",el->e);
			fclose(file);
}

int main(int argc,char** argv){
	edgelist* el;
	graph* g;
	unsigned char k=atoi(argv[2]);
	unsigned long long nck;
	unsigned test = 0;
	if (argc > 4) {
		test = 1;
	}
	omp_set_num_threads(atoi(argv[1]));

	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;

	if (!test) {
		printf("Reading edgelist from file %s\n",argv[3]);
	}

	el=readedgelist(argv[3]);

	if (!test) {
		printf("Number of nodes = %u\n",el->n);
		printf("Number of edges = %u\n",el->e);
	}
	unsigned number_of_nodes = el->n;


	t2=time(NULL);
	if (!test) {
		printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	}


	if (!test) {
		printf("Building the graph structure\n");
	}
	ord_core(el);
	relabel(el);
	g=mkgraph(el);

	if (!test) {
		printf("Number of nodes (degree > 0) = %u\n",g->n);
	}

	t2=time(NULL);
	if (!test) {
		printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	}
	t1=t2;

	unsigned i,n_m,e_m;
	bool *rm=(bool*)calloc(g->n,sizeof(bool));
	allocglobal(g,k,number_of_nodes);//allocataing global variables
	nck=kclique_main(k, g, el->node_map);

	if (!test){
		printf("Number of %u-cliques: %llu\n",k,nck);
	}
	unsigned long long r=0,r2=0;
	free_graph(g);

	t2=time(NULL);
	t1=t2;

	//iterate through clique_matrix to print
	for (int p = 0; p < ck_m->vector_length; p++) {//vector iterator
		for (std::pair<unsigned, int64_t> element : clique_matrix->clique_mat[p]) {
				printf("%u %u %ld\n", p, element.first, element.second);
		}
	}

	if (!test) {
		printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));
	}

	return 0;
}
