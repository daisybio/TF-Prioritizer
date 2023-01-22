import { ComponentFixture, TestBed } from '@angular/core/testing';

import { CoOccurrenceComponent } from './co-occurrence.component';

describe('CoOccurrenceComponent', () => {
  let component: CoOccurrenceComponent;
  let fixture: ComponentFixture<CoOccurrenceComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ CoOccurrenceComponent ]
    })
    .compileComponents();

    fixture = TestBed.createComponent(CoOccurrenceComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
