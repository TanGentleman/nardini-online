import React from 'react';
import { useReactTable, getCoreRowModel, getSortedRowModel, flexRender } from '@tanstack/react-table';
import type { ColumnDef } from '@tanstack/react-table';
import type { SequenceData } from '../types/schemas';
import { StatusBadge } from './StatusBadge';

interface SequencesTableProps {
  sequences: SequenceData[];
}

export const SequencesTable: React.FC<SequencesTableProps> = ({ sequences }) => {
  const columns: ColumnDef<SequenceData>[] = [
    {
      accessorKey: 'sequence_id',
      header: 'Sequence ID',
      cell: ({ getValue }) => (
        <span className="font-mono text-sm font-medium text-gray-900">
          {getValue() as string}
        </span>
      ),
    },
    {
      accessorKey: 'status',
      header: 'Status',
      cell: ({ getValue }) => (
        <StatusBadge status={getValue() as any} type="sequence" />
      ),
    },
    {
      accessorKey: 'seq_uuid',
      header: 'UUID',
      cell: ({ getValue }) => (
        <span className="font-mono text-xs text-gray-500">
          {(getValue() as string)?.slice(0, 8)}...
        </span>
      ),
    },
    {
      accessorKey: 'start_time',
      header: 'Started',
      cell: ({ getValue }) => {
        const timestamp = getValue() as number;
        return (
          <span className="text-sm text-gray-600">
            {timestamp ? new Date(timestamp * 1000).toLocaleString() : '—'}
          </span>
        );
      },
    },
    {
      accessorKey: 'end_time',
      header: 'Completed',
      cell: ({ getValue }) => {
        const timestamp = getValue() as number;
        return (
          <span className="text-sm text-gray-600">
            {timestamp ? new Date(timestamp * 1000).toLocaleString() : '—'}
          </span>
        );
      },
    },
    {
      accessorKey: 'job_id',
      header: 'Job ID',
      cell: ({ getValue }) => (
        <span className="font-mono text-xs text-gray-500">
          {(getValue() as string) || '—'}
        </span>
      ),
    },
    {
      accessorKey: 'zip_path',
      header: 'Download',
      cell: ({ getValue }) => {
        const zipPath = getValue() as string;
        return zipPath ? (
          <a
            href={`/api/download/${zipPath}`}
            className="text-blue-600 hover:text-blue-800 transition-colors"
            download
          >
            📁
          </a>
        ) : (
          <span className="text-gray-400">—</span>
        );
      },
    },
  ];

  const table = useReactTable({
    data: sequences,
    columns,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
  });

  return (
    <div className="bg-white rounded-lg border border-gray-200 overflow-hidden">
      <div className="px-6 py-4 border-b border-gray-200">
        <h3 className="text-lg font-semibold text-gray-900">Sequences</h3>
      </div>

      <div className="overflow-x-auto">
        <table className="min-w-full divide-y divide-gray-200">
          <thead className="bg-gray-50">
            {table.getHeaderGroups().map(headerGroup => (
              <tr key={headerGroup.id}>
                {headerGroup.headers.map(header => (
                  <th
                    key={header.id}
                    className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider cursor-pointer hover:bg-gray-100"
                    onClick={header.column.getToggleSortingHandler()}
                  >
                    {header.isPlaceholder ? null : (
                      <div className="flex items-center space-x-1">
                        <span>{flexRender(header.column.columnDef.header, header.getContext())}</span>
                        <span className="text-gray-400">
                          {header.column.getIsSorted() === 'asc' ? '↑' :
                           header.column.getIsSorted() === 'desc' ? '↓' : '↕'}
                        </span>
                      </div>
                    )}
                  </th>
                ))}
              </tr>
            ))}
          </thead>
          <tbody className="bg-white divide-y divide-gray-200">
            {table.getRowModel().rows.length === 0 ? (
              <tr>
                <td colSpan={columns.length} className="px-6 py-12 text-center text-gray-500">
                  No sequences found
                </td>
              </tr>
            ) : (
              table.getRowModel().rows.map(row => (
                <tr key={row.id} className="hover:bg-gray-50 transition-colors">
                  {row.getVisibleCells().map(cell => (
                    <td key={cell.id} className="px-6 py-4 whitespace-nowrap">
                      {flexRender(cell.column.columnDef.cell, cell.getContext())}
                    </td>
                  ))}
                </tr>
              ))
            )}
          </tbody>
        </table>
      </div>
    </div>
  );
};