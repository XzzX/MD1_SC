#include    <iostream>

#include	"CellSubdivision.h"
#include	"particles.hpp"

/*
@param width width of the box
@param height height of the box
@param rc cut-off distance
*/
CellSubdivision::CellSubdivision(const double width, const double height, const double rc):
	m_origin(0.0,0.0,0.0), m_boxHeight(height), m_boxWidth(width), m_numCellsX(floor(width/rc)), m_numCellsY(floor(height/rc)){

	m_cellWidth = (width/static_cast<double>(m_numCellsX));
	m_cellHeight = (height/static_cast<double>(m_numCellsY));
	//create cell array
	m_cells.resize(m_numCellsX*m_numCellsY);

}

/*
@param pos position of particle
@param id id of particle
@return cell id
*/
unsigned int	CellSubdivision::InsertParticle(const Vector& pos, const unsigned int id){
	Vector a(pos - m_origin);
	unsigned int x(floor(a[0]/m_cellWidth));
	unsigned int y(floor(a[1]/m_cellHeight));
	unsigned int cellId(GetCell(x,y));
	m_cells[cellId].insert(id);
	return cellId;
}

/*
@param cellId id of the cell the particle is in
@param id id of the particle
*/
void	CellSubdivision::DeleteParticle(const unsigned int cellId, const unsigned int id){
	m_cells[cellId].erase(id);
}

/*
@param cellId id of the cell the particle is in
@param id id of the particle
@return list iterator of a list of nearest neighbours
*/
list<unsigned int>::iterator CellSubdivision::GetNeighbours(const unsigned int cellId, const unsigned int id){
    m_neighbours.clear();

    int x,y;
    int a,b;
    GetCoordinates(cellId, x, y);

    set<unsigned int>::iterator it;

    for (int i = x-1; i<x+2; i++)
        for (int j = y-1; j<y+2; j++){
            a = i;
            b = j;
            NormaliseCoordinates(a,b);
            for (it=m_cells[GetCell(a,b)].begin(); it!=m_cells[GetCell(a,b)].end(); it++)
                if (*it!=id) m_neighbours.push_back(*it);
        }

    return m_neighbours.begin();
}
