import React from 'react';
import type { SequenceStatus, RunStatus } from '../types/schemas';

interface StatusBadgeProps {
  status: SequenceStatus | RunStatus;
  type: 'sequence' | 'run';
}

const statusConfig = {
  sequence: {
    complete: { bg: 'bg-green-100', text: 'text-green-800', icon: '✅' },
    pending: { bg: 'bg-yellow-100', text: 'text-yellow-800', icon: '⏳' },
    cached: { bg: 'bg-blue-100', text: 'text-blue-800', icon: '💾' },
    pending_external: { bg: 'bg-orange-100', text: 'text-orange-800', icon: '🔄' }
  },
  run: {
    complete: { bg: 'bg-green-100', text: 'text-green-800', icon: '✅' },
    pending: { bg: 'bg-yellow-100', text: 'text-yellow-800', icon: '⏳' }
  }
};

export const StatusBadge: React.FC<StatusBadgeProps> = ({ status, type }) => {
  const config = statusConfig[type][status as keyof typeof statusConfig[typeof type]];

  return (
    <span className={`inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium ${config.bg} ${config.text}`}>
      <span className="mr-1">{config.icon}</span>
      {status}
    </span>
  );
};