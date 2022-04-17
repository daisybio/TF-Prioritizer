import { TestBed } from '@angular/core/testing';

import { MapToTableService } from './map-to-table.service';

describe('MapToTableService', () => {
  let service: MapToTableService;

  beforeEach(() => {
    TestBed.configureTestingModule({});
    service = TestBed.inject(MapToTableService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
