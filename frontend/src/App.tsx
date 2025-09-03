import React, { useState, useEffect } from 'react';
import type { RunData } from './types/schemas';
import { mockRuns } from './data/mockData';
import { SequencesTable } from './components/SequencesTable';
import { StatusBadge } from './components/StatusBadge';
import './App.css';

// API configuration
const API_BASE_URL = 'http://localhost:8000';

const App: React.FC = () => {
  const [runs, setRuns] = useState<RunData[]>([]);
  const [selectedRun, setSelectedRun] = useState<RunData | null>(null);
  const [loading, setLoading] = useState(false);
  const [importing, setImporting] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const fetchRuns = async () => {
    try {
      setLoading(true);
      setError(null);
      const response = await fetch(`${API_BASE_URL}/runs`);
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
      const data = await response.json() as { runs: RunData[] };
      setRuns(data.runs || []);
    } catch (err) {
      console.warn('Using mock data due to API error:', err);
      setRuns(mockRuns);
    } finally {
      setLoading(false);
    }
  };

  const importData = async () => {
    try {
      setImporting(true);
      setError(null);
      await new Promise(resolve => setTimeout(resolve, 2000));
      await fetchRuns();
    } catch (err) {
      setError('Failed to import data');
    } finally {
      setImporting(false);
    }
  };

  useEffect(() => {
    fetchRuns();
  }, []);

  if (loading) {
    return (
      <div className="min-h-screen bg-gray-50 flex items-center justify-center">
        <div className="text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto"></div>
          <p className="mt-4 text-gray-600">Loading runs...</p>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-50">
      <div className="max-w-7xl mx-auto px-4 py-8">
        {/* Header */}
        <div className="flex justify-between items-center mb-8">
          <div>
            <h1 className="text-3xl font-bold text-gray-900">Nardini Sequence Analysis</h1>
            <p className="text-gray-600 mt-1">Monitor and download sequence processing results</p>
          </div>
          <button
            onClick={importData}
            disabled={importing}
            className="inline-flex items-center px-4 py-2 border border-transparent text-sm font-medium rounded-md text-white bg-blue-600 hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500 disabled:opacity-50"
          >
            {importing ? (
              <>
                <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-white mr-2"></div>
                Importing...
              </>
            ) : (
              <>
                <span className="mr-2">📥</span>
                Import Data
              </>
            )}
          </button>
        </div>

        {/* Error Alert */}
        {error && (
          <div className="mb-6 bg-red-50 border border-red-200 rounded-md p-4">
            <div className="flex">
              <div className="flex-shrink-0">
                <span className="text-red-400">⚠️</span>
              </div>
              <div className="ml-3">
                <p className="text-sm text-red-700">{error}</p>
              </div>
            </div>
          </div>
        )}

        {/* Runs List */}
        <div className="bg-white rounded-lg border border-gray-200 overflow-hidden mb-8">
          <div className="px-6 py-4 border-b border-gray-200">
            <h2 className="text-xl font-semibold text-gray-900">Analysis Runs</h2>
            <p className="text-gray-600 text-sm">Available sequence analysis runs</p>
          </div>

          <div className="divide-y divide-gray-200">
            {runs.length === 0 ? (
              <div className="px-6 py-12 text-center text-gray-500">
                <span className="text-4xl mb-4 block">📋</span>
                <p>No runs available. Click "Import Data" to load runs.</p>
              </div>
            ) : (
              runs.map(run => (
                <div key={run.submitted_at} className="px-6 py-4 hover:bg-gray-50 transition-colors">
                  <div className="flex justify-between items-start">
                    <div className="flex-1">
                      <div className="flex items-center space-x-3 mb-2">
                        <h3 className="text-lg font-medium text-gray-900">
                          {run.fasta_filename}
                        </h3>
                        <StatusBadge status={run.status} type="run" />
                      </div>
                      <div className="text-sm text-gray-600 space-y-1">
                        <p>Sequences: {run.total_sequences} total, {run.cached_sequences} cached</p>
                        <p>Submitted: {new Date(run.submitted_at * 1000).toLocaleString()}</p>
                        {run.completed_at && (
                          <p>Completed: {new Date(run.completed_at * 1000).toLocaleString()}</p>
                        )}
                      </div>
                    </div>
                    <div className="flex space-x-3 ml-4">
                      <button
                        onClick={() => setSelectedRun(run)}
                        className="inline-flex items-center px-3 py-2 border border-gray-300 shadow-sm text-sm leading-4 font-medium rounded-md text-gray-700 bg-white hover:bg-gray-50 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500"
                      >
                        <span className="mr-2">👁️</span>
                        View Sequences
                      </button>
                      {run.merged_zip_filename && (
                        <a
                          href={`/api/download/${run.merged_zip_filename}`}
                          className="inline-flex items-center px-3 py-2 border border-transparent text-sm leading-4 font-medium rounded-md text-white bg-green-600 hover:bg-green-700 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-green-500"
                          download
                        >
                          <span className="mr-2">📁</span>
                          Download All
                        </a>
                      )}
                    </div>
                  </div>
                </div>
              ))
            )}
          </div>
        </div>

        {/* Sequences Table */}
        {selectedRun && (
          <div className="space-y-6">
            <div className="flex justify-between items-center">
              <div>
                <h2 className="text-2xl font-bold text-gray-900">Sequences</h2>
                <p className="text-gray-600">{selectedRun.fasta_filename}</p>
              </div>
              <button
                onClick={() => setSelectedRun(null)}
                className="inline-flex items-center px-3 py-2 border border-gray-300 shadow-sm text-sm leading-4 font-medium rounded-md text-gray-700 bg-white hover:bg-gray-50 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-gray-500"
              >
                <span className="mr-2">✕</span>
                Close
              </button>
            </div>

            <SequencesTable sequences={Object.values(selectedRun.sequences)} />
          </div>
        )}
      </div>
    </div>
  );
};

export default App;