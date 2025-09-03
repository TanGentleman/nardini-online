# Nardini Online Frontend

🖥️ **React frontend for Nardini Online protein analysis**

A modern React application that provides a user-friendly interface for monitoring and downloading protein sequence analysis results from the Nardini Online backend.

## 🚀 Features

- **Real-time Monitoring**: View analysis runs and their current status
- **Sequence Visualization**: Interactive table displaying processed sequences
- **Data Import**: Import new analysis runs from the backend
- **File Downloads**: Download merged analysis results as ZIP files
- **Responsive Design**: Clean, modern UI built with Tailwind CSS

## 🛠️ Technology Stack

- **Runtime**: [Bun](https://bun.com) - Fast all-in-one JavaScript runtime
- **Framework**: React 19 with TypeScript
- **Styling**: Tailwind CSS 4
- **Tables**: TanStack React Table
- **Build Tool**: Vite (via Bun)

## 📁 Project Structure

```
frontend/
├── src/
│   ├── components/
│   │   ├── SequencesTable.tsx    # Interactive sequences table
│   │   └── StatusBadge.tsx       # Status indicator component
│   ├── types/
│   │   └── schemas.ts            # TypeScript type definitions
│   ├── App.tsx                   # Main application component
│   └── index.tsx                 # Application entry point
├── dist/                         # Built application files
├── package.json
├── tailwind.config.js
├── postcss.config.js
└── tsconfig.json
```

## 🚀 Getting Started

### Prerequisites

- [Bun](https://bun.sh) runtime installed
- Backend API running on `http://localhost:8000`

### Installation

```bash
# Install dependencies
bun install
```

### Development

```bash
# Start development server with hot reload
bun run dev
```

The application will be available at `http://localhost:8080` (or the next available port).

### Build

```bash
# Build for production
bun run build
```

### Serve

```bash
# Serve built files locally
bun run serve
```

## 🔧 Configuration

The frontend connects to the backend API at `http://localhost:8000` by default. To change this:

1. Open `src/App.tsx`
2. Update the `API_BASE_URL` constant:

```typescript
const API_BASE_URL = 'http://your-backend-url:port';
```

## 📊 Usage

1. **View Analysis Runs**: The main dashboard displays all available analysis runs
2. **Monitor Status**: Status badges show the current state of each run
3. **Import Data**: Click "Import Data" to refresh runs from the backend
4. **View Sequences**: Click "View Sequences" to see detailed sequence data in a table
5. **Download Results**: Download merged analysis results as ZIP files

## 🔗 Integration

This frontend is designed to work with the Nardini Online backend API. Make sure the backend is running and accessible before starting the frontend.

For backend setup instructions, see the main project [README](../README.md) and [GUIDE](../GUIDE.md).

## 📝 Development Notes

- Built with React 19 and modern hooks
- Uses TypeScript for type safety
- Styled with Tailwind CSS for consistent design
- Responsive design works on desktop and mobile devices
- Hot reload enabled during development for fast iteration
