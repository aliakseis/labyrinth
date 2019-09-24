#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>
#include <queue>

#include <assert.h>

using std::cout;
using std::cin;
using std::cerr;
using std::ifstream;
using std::exception;
using std::string;
using std::vector;



struct Plex     // Similar to MFC CPlex
	// warning variable length structure
{
	Plex* pNext;
	int dwReserved[1];    // align on 8 byte boundary

	void* data() { return this + 1; }

	// like 'calloc' but no zero fill
	// may throw memory exceptions
	static Plex* Create(Plex*& pHead, size_t nMax, size_t cbElement)
	{
		assert(nMax > 0 && cbElement > 0);
		Plex* p = (Plex*) operator new(sizeof(Plex) + nMax * cbElement);
		// may throw exception
		p->pNext = pHead;
		pHead = p;  // change head (adds in reverse order for simplicity)
		return p;
	}

	void FreeDataChain()       // free this one and links
	{
		Plex* p = this;
		while (p != 0)
		{
			void* bytes = p;
			p = p->pNext;
			operator delete(bytes);
		}
	}
};

typedef unsigned int UINT;

class FixedAlloc    // Similar to MFC CFixedAlloc
{
	// Constructors
public:
	FixedAlloc(UINT nAllocSize, UINT nBlockSize = 64);

	// Attributes
	UINT GetAllocSize() { return m_nAllocSize; }

	// Operations
public:
	void* Alloc();  // return a chunk of memory of nAllocSize
	void Free(void* p); // free chunk of memory returned from Alloc
	void FreeAll(); // free everything allocated from this allocator

	// Implementation
public:
	~FixedAlloc();

protected:
	struct CNode
	{
		CNode* pNext;	// only valid when in free list
	};

	UINT m_nAllocSize;	// size of each block from Alloc
	UINT m_nBlockSize;	// number of blocks to get at a time
	Plex* m_pBlocks;	// linked list of blocks (is nBlocks*nAllocSize)
	CNode* m_pNodeFree;	// first free node (0 if no free nodes)
};

FixedAlloc::FixedAlloc(UINT nAllocSize, UINT nBlockSize)
{
	assert(nAllocSize >= sizeof(CNode));
	assert(nBlockSize > 1);

	if (nAllocSize < sizeof(CNode))
		nAllocSize = sizeof(CNode);
	if (nBlockSize <= 1)
		nBlockSize = 64;

	m_nAllocSize = nAllocSize;
	m_nBlockSize = nBlockSize;
	m_pNodeFree = 0;
	m_pBlocks = 0;
}

FixedAlloc::~FixedAlloc()
{
	FreeAll();
}

void FixedAlloc::FreeAll()
{
	m_pBlocks->FreeDataChain();
	m_pBlocks = 0;
	m_pNodeFree = 0;
}

void* FixedAlloc::Alloc()
{
	if (m_pNodeFree == 0)
	{
		// add another block
		Plex* pNewBlock = Plex::Create(m_pBlocks, m_nBlockSize, m_nAllocSize);

		// chain them into free list
		// free in reverse order to make it easier to debug
		char* pData = (char*)pNewBlock->data() + m_nAllocSize * (m_nBlockSize - 1);
		for (int i = m_nBlockSize; i > 0; i--, pData -= m_nAllocSize)
		{
			CNode* pNode = (CNode*)pData;
			pNode->pNext = m_pNodeFree;
			m_pNodeFree = pNode;
		}
	}
	assert(m_pNodeFree != 0);  // we must have something

	// remove the first available node from the free list
	void* pNode = m_pNodeFree;
	m_pNodeFree = m_pNodeFree->pNext;
	return pNode;
}

void FixedAlloc::Free(void* p)
{
	if (p != 0)
	{
		// simply return the node to the free list
		CNode* pNode = (CNode*)p;
		pNode->pNext = m_pNodeFree;
		m_pNodeFree = pNode;
	}
}

//////////////////////////////////////////////////////////////////////////////


enum { HASH_BITS = 12, HASH_SIZE = 1 << HASH_BITS };


int end_x, end_y;

struct VertexKey
{
	int x, y;

	int hash() const
	{
		return (unsigned int)(((x << 1) ^ y) * 2654435769U) >> (32 - HASH_BITS);
	}

	// Manhattan
	unsigned int distance() const { return abs(x - end_x) + abs(y - end_y); }

};

inline bool operator == (const VertexKey& left, const VertexKey& right)
{
	return left.x == right.x && left.y == right.y;
}


struct Vertex
{
	Vertex* next;
	VertexKey key;
	union
	{
		std::size_t pathLength; // must be odd
		const Vertex* parent;
	};

	bool empty() const { return key.x == 0 && key.y == 0; }

	bool visited() const { return (pathLength & 1u) == 0; }
	// Manhattan
	//unsigned int estimate() { return pathLength + abs(key.x - end_x) + abs(key.y - end_y); }


	void* operator new(size_t size)
	{
		assert(size == s_alloc.GetAllocSize());
		return s_alloc.Alloc();
	}
	void* operator new(size_t, void* p){ return p; }
	void operator delete(void* p) { s_alloc.Free(p); }

protected:
	static FixedAlloc s_alloc;
};

FixedAlloc Vertex::s_alloc(sizeof(Vertex), 64);

Vertex hashTable[HASH_SIZE];


class Priority
{
	friend bool operator < (const Priority& left, const Priority& right);
public:
	Vertex* vertex;
	Vertex* parent;
	unsigned int weight() { return weight_finishing >> 1; }
	bool finishing() { return !(weight_finishing & 1u); }
	Priority(Vertex* vertex_, Vertex* parent_, unsigned int weight_, bool finishing_)
		: vertex(vertex_), parent(parent_), weight_finishing((weight_ << 1) + !finishing_) {}
private:
	unsigned int weight_finishing;
};

inline bool operator < (const Priority& left, const Priority& right)
{
	return left.weight_finishing > right.weight_finishing;
}



Vertex* findItem(const VertexKey& newKey, int idx)
{
	Vertex* vertex = hashTable + idx;
	if (!vertex->empty())
	{
		do
		{
			if (vertex->key == newKey)
			{
				return vertex;
			}
			vertex = vertex->next;
		} while (vertex != 0);
	}

	return 0;
}


typedef std::priority_queue<Priority> PriorityQueue;



vector<string> labyrinth;


void reportSolution(const Vertex* vertex, const Vertex* parent)
{
	for (;;)
	{
		labyrinth[vertex->key.y][vertex->key.x] = '+';
		vertex = parent;
		if (!vertex)
			break;
		parent = vertex->parent;
	}

	for (const auto& s : labyrinth)
	{
		cout << s << '\n';
	}
}

unsigned int minFinishing = 0xFFFFFFFF;

const int offsets[] = { 0, 1, 0, -1, 0 };

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		cerr << "Usage: test input_file\n";
		return 1;
	}

	try
	{
		ifstream inFile(argv[1]);
		string str;
		while (getline(inFile, str))
		{
			if (str.empty())
				break;

			labyrinth.push_back(str);
		}

		if (labyrinth.empty())
			return 0;


		PriorityQueue queue;

		const int start_x = labyrinth[0].find(' ');
		end_x = labyrinth[end_x].find(' ');
		end_y = labyrinth.size() - 1;

		labyrinth[0][start_x] = '+';

		if (labyrinth.size() == 1)
		{
			cout << labyrinth[0] << '\n';
			return 0;
		}

		VertexKey initialKey = { start_x, 1 };
		int initialIdx = initialKey.hash();
		Vertex* v = hashTable + initialIdx;
		v->next = 0;
		v->key = initialKey;
		v->pathLength = 1; // must be odd
		queue.push(Priority(v, 0, initialKey.distance() + 1, false));

		do
		{
			Priority priority = queue.top();
			if (priority.finishing())
			{
				reportSolution(priority.vertex, priority.parent);
				break;
			}

			queue.pop();

			if (priority.vertex->visited())
			{
				continue;
			}


			for (int direction = 4; --direction >= 0; )

			{
				const VertexKey newKey = { 
					priority.vertex->key.x + offsets[direction], 
					priority.vertex->key.y + offsets[direction + 1]
				};

				if (labyrinth[newKey.y][newKey.x] != ' ')
				{
					continue;
				}

				const int idx = newKey.hash();
				Vertex* v = findItem(newKey, idx);
				if (v && v->visited())
				{
					continue;
				}

				if (v && v->pathLength <= priority.vertex->pathLength + 2)
				{
					continue;
				}

				const bool finishing = newKey.x == end_x && newKey.y == end_y;

				const unsigned int estimate = newKey.distance() + (priority.vertex->pathLength + 2) / 2;
				if (estimate >= minFinishing)
				{
					continue;
				}

				if (finishing)
				{
					if (minFinishing > estimate)
					{
						minFinishing = estimate;
					}
				}


				if (!v)
				{
					if (hashTable[idx].empty())
					{
						v = hashTable + idx;
					}
					else
					{
						v = new Vertex();
						v->next = hashTable[idx].next;
						hashTable[idx].next = v;
					}
					v->key = newKey;
				}

				v->pathLength = priority.vertex->pathLength + 2; // must be odd

				queue.push(Priority(v, priority.vertex, estimate, finishing));
			}

			// visited
			priority.vertex->parent = priority.parent;
		} while (!queue.empty());


	}
	catch (exception& e)
	{
		cerr << e.what();
		return 1;
	}

	return 0;
}
